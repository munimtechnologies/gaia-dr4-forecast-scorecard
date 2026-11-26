# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
# Made with AI agents working in MV2's name.
"""Tests for the Gaia DR4 forecast scorecard (standard library only, no network)."""
import json
import math
import os
import tempfile
import unittest

import scorecard as sc
import job

HERE = os.path.dirname(os.path.abspath(__file__))


class Rules(unittest.TestCase):
    def test_parse_periods(self):
        self.assertEqual(sc.parse_periods("3054.9"), [(3054.9, 3054.9)])
        self.assertEqual(sc.parse_periods("1700;3000"), [(1700.0, 1700.0), (3000.0, 3000.0)])
        self.assertEqual(sc.parse_periods("200-300"), [(200.0, 300.0)])
        self.assertEqual(sc.parse_periods(""), [])

    def test_period_match(self):
        self.assertTrue(sc.period_match(2500, "3054.9"))
        self.assertFalse(sc.period_match(2400, "3054.9"))
        self.assertTrue(sc.period_match(3600, "3054.9"))
        self.assertFalse(sc.period_match(3700, "3054.9"))
        self.assertTrue(sc.period_match(2900, "1700;3000"))
        self.assertTrue(sc.period_match(350, "200-300"))
        self.assertFalse(sc.period_match(370, "200-300"))
        self.assertFalse(sc.period_match(None, "256"))

    def test_classify(self):
        self.assertEqual(sc.classify_solution("OrbitalTargetedSearchValidated", "gaiadr4.nss_two_body_orbit"), "orbit")
        self.assertEqual(sc.classify_solution("AstroSpectroSB1", "gaiadr4.nss_two_body_orbit"), "orbit")
        self.assertEqual(sc.classify_solution("Acceleration9", "gaiadr4.nss_acceleration_astro"), "accel")
        self.assertEqual(sc.classify_solution("SB1", "gaiadr4.nss_two_body_orbit"), "other")
        self.assertEqual(sc.classify_solution("FirstDegreeTrendSB1", "gaiadr4.nss_non_linear_spectro"), "other")
        self.assertEqual(sc.classify_solution("EclipsingBinary", "gaiadr4.nss_two_body_orbit"), "other")
        self.assertEqual(sc.classify_solution("", "gaiadr4.exoplanet_candidates", True), "orbit")
        self.assertEqual(sc.classify_solution("", "gaiadr4.nss_acceleration_astro"), "accel")
        self.assertEqual(sc.classify_solution("", "gaiadr4.nss_non_linear_spectro"), "other")

    def test_host_status_and_labels(self):
        self.assertEqual(sc.host_status([], "100"), "NOTHING")
        orb = {"kind": "orbit", "period": 100.0}
        self.assertEqual(sc.host_status([orb], "3054.9"), "ORBIT")
        self.assertEqual(sc.host_status([orb, {"kind": "orbit", "period": 3000.0}], "3054.9"), "ORBIT_MATCH")
        self.assertEqual(sc.host_status([{"kind": "accel", "period": None}], "3054.9"), "ACCEL")
        self.assertEqual(sc.host_status([{"kind": "other", "period": 3054.9}], "3054.9"), "NOTHING")
        self.assertEqual(sc.label_detectable("ORBIT_MATCH"), "RIGHT")
        self.assertEqual(sc.label_detectable("ORBIT"), "PARTLY")
        self.assertEqual(sc.label_detectable("ACCEL"), "PARTLY")
        self.assertEqual(sc.label_detectable("NOTHING"), "WRONG")
        self.assertEqual(sc.label_not_detectable("ORBIT_MATCH"), "WRONG")
        for s in ("ORBIT", "ACCEL", "NOTHING"):
            self.assertEqual(sc.label_not_detectable(s), "RIGHT")
        self.assertEqual(sc.label_for(sc.DETECTABLE, None), "UNGRADED_NO_COUNTERPART")

    def test_probabilistic(self):
        self.assertAlmostEqual(sc.brier(0.435, 1), 0.319225)
        self.assertAlmostEqual(sc.brier(0.435, 0), 0.189225)
        self.assertEqual(sc.w1_outcome("1.4"), 1)
        self.assertEqual(sc.w1_outcome(1.15), 0)
        self.assertIsNone(sc.w1_outcome(""))
        y, dra, ddec = sc.w2_outcome(10.2, -5.01, 10.0, -5.0)
        self.assertEqual(y, 1)
        self.assertAlmostEqual(dra, 0.2)
        y, _, _ = sc.w2_outcome(10.05, -5.05, 10.0, -5.0)
        self.assertEqual(y, 0)
        y, _, _ = sc.w2_outcome(10.0, -5.1, 10.0, -5.0)
        self.assertEqual(y, 1)
        self.assertEqual(sc.w2_outcome(None, 1, 1, 1)[0], None)

    def test_thiele_innes_round_trip(self):
        for inc in (11.0, 60.0, 89.0, 120.0, 170.0):
            for w, W in ((30.0, 50.0), (200.0, 310.0)):
                A, B, F, G = sc.thiele_innes_forward(0.7, inc, w, W)
                self.assertAlmostEqual(sc.inclination_from_thiele_innes(A, B, F, G), inc, places=6)
        self.assertIsNone(sc.inclination_from_thiele_innes(None, 1, 1, 1))

    def test_b3(self):
        self.assertEqual(sc.b3_label("ORBIT_MATCH", 169.0), "BH17_SUPPORTED")
        self.assertEqual(sc.b3_label("ORBIT_MATCH", 7.0), "BH17_SUPPORTED")
        self.assertEqual(sc.b3_label("ORBIT_MATCH", 45.0), "NOT_BH17")
        self.assertEqual(sc.b3_label("ORBIT_MATCH", 150.0), "NOT_BH17")
        self.assertEqual(sc.b3_label("ORBIT_MATCH", 20.0), "INTERMEDIATE")
        self.assertEqual(sc.b3_label("ORBIT_MATCH", None), "MATCH_NO_INCLINATION")
        self.assertEqual(sc.b3_label("NOTHING", None), "NO_ORBIT")
        self.assertEqual(sc.b3_label(None, None), "NO_ORBIT")
        self.assertEqual(sc.b3_label("ORBIT", 10.0), "NO_MATCHING_ORBIT")

    def test_l1(self):
        self.assertEqual(sc.l1_label(3300), "RIGHT")
        self.assertEqual(sc.l1_label(11700), "RIGHT")
        self.assertEqual(sc.l1_label(3299), "LOW")
        self.assertEqual(sc.l1_label(11701), "HIGH")
        self.assertEqual(sc.l1_label(None), "UNGRADABLE")

    def test_crossmatch(self):
        dr3 = {"source_id": "42", "ra": 10.0, "dec": 20.0, "pmra": 0.0, "pmdec": 1000.0, "phot_g_mean_mag": 9.0}
        moved = {"source_id": "42", "ra": 10.0, "dec": 20.0 + 1.5 / 3600, "phot_g_mean_mag": 9.1}
        old = {"source_id": "43", "ra": 10.0, "dec": 20.0, "phot_g_mean_mag": 9.0}
        bright = {"source_id": "44", "ra": 10.0, "dec": 20.0 + 1.5 / 3600, "phot_g_mean_mag": 8.4}
        far = {"source_id": "45", "ra": 10.0, "dec": 20.0 + 4.0 / 3600, "phot_g_mean_mag": 9.0}
        best, info = sc.crossmatch_pick(dr3, [old, bright, far, moved], 2017.5, 2016.0)
        self.assertEqual(best["source_id"], "42")
        self.assertLess(info["sep_arcsec"], 0.01)
        self.assertEqual(info["n_candidates_ok"], 2)
        self.assertTrue(info["same_source_id"])
        best, info = sc.crossmatch_pick(dr3, [far, bright], 2017.5)
        self.assertIsNone(best)

    def test_helpers(self):
        self.assertEqual(sc.id_str("Gaia DR3 5287961677550103424"), "5287961677550103424")
        self.assertEqual(sc.id_str(None), "")
        self.assertEqual(sc.id_str("1.5"), "")
        rows = sc.read_csv_text("# note\nA,b\n1,x\n\n2,y\n")
        self.assertEqual(rows, [{"a": "1", "b": "x"}, {"a": "2", "b": "y"}])
        self.assertIsNone(sc.to_float("null"))
        self.assertEqual(sc.to_float("2.5"), 2.5)


class FrozenTargets(unittest.TestCase):
    def test_targets_match_preregistration(self):
        t = job.load_targets(os.path.join(HERE, "targets.csv"))
        self.assertEqual(len(t), 24)
        self.assertEqual(len({x["id"] for x in t}), 24)
        groups = {}
        for x in t:
            groups[x["group"]] = groups.get(x["group"], 0) + 1
        self.assertEqual(groups, {"detectable": 9, "not_detectable": 11, "w1": 1, "w2": 1, "b3": 1, "l1": 1})
        det = sorted(x["id"] for x in t if x["group"] == "detectable")
        self.assertEqual(det, sorted(["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P18", "B1"]))
        for x in t:
            if x["group"] != "l1":
                self.assertTrue(sc.id_str(x["dr3_source_id"]), x["id"])
                self.assertTrue(-90 <= float(x["dec"]) <= 90)
            if x["group"] in ("detectable", "not_detectable", "b3"):
                self.assertTrue(sc.parse_periods(x["period_spec"]), x["id"])
        p = {x["id"]: x for x in t}
        self.assertEqual(p["W1"]["p"], "0.435")
        self.assertEqual(p["W2"]["p"], "0.948")
        self.assertEqual(p["B1"]["period_spec"], "1823.5")
        self.assertEqual(p["B2"]["period_spec"], "228.8")

    def test_io_module_imports(self):
        import io_dr4score
        for name in ("preflight", "list_tables", "list_columns", "dr4_cone", "rows_by_ids", "count",
                     "dr3_rows", "simbad_dr3"):
            self.assertTrue(callable(getattr(io_dr4score, name)))


# ---------------------------------------------------------------- end-to-end with a fake archive
IDS = {k: "10000000000000000%02d" % i for i, k in enumerate("ABCDEFG", 1)}


def _tables():
    A, B, C, D, E, F = (IDS[k] for k in "ABCDEF")
    TI = sc.thiele_innes_forward(0.9, 11.0, 30.0, 50.0)
    return {
        "gaiadr4.nss_two_body_orbit": [
            {"source_id": A, "nss_solution_type": "OrbitalTargetedSearch", "period": "2900"},
            {"source_id": B, "nss_solution_type": "SB1", "period": "10"},
            {"source_id": C, "nss_solution_type": "Orbital", "period": "430"},
            {"source_id": F, "nss_solution_type": "Orbital", "period": "255.6",
             "a_thiele_innes": str(TI[0]), "b_thiele_innes": str(TI[1]),
             "f_thiele_innes": str(TI[2]), "g_thiele_innes": str(TI[3])},
        ],
        "gaiadr4.nss_acceleration_astro": [{"source_id": D, "nss_solution_type": "Acceleration7"}],
        "gaiadr4.exoplanet_candidates": [],
        "gaiadr4.dr3_neighbourhood": [{"dr3_source_id": IDS[k], "dr4_source_id": IDS[k]} for k in "ABCDEF"],
    }


COLS = {
    "gaiadr4.nss_two_body_orbit": ["source_id", "nss_solution_type", "period", "a_thiele_innes",
                                   "b_thiele_innes", "f_thiele_innes", "g_thiele_innes"],
    "gaiadr4.nss_acceleration_astro": ["source_id", "nss_solution_type"],
    "gaiadr4.exoplanet_candidates": ["source_id", "period", "mass"],
    "gaiadr4.dr3_neighbourhood": ["dr3_source_id", "dr4_source_id"],
}


class FakeIO:
    def __init__(self, ready=True):
        self.ready = ready
        self.tables = _tables()
        self.pos = {IDS[k]: (10.0 * i, -20.0 + i) for i, k in enumerate("ABCDEFG", 1)}

    def preflight(self, schema):
        return self.ready, "fake"

    def list_tables(self, schema):
        names = ["gaiadr4.gaia_source"] + list(self.tables)
        return [{"schema_name": "gaiadr4", "table_name": n, "description": ""} for n in names]

    def list_columns(self, table):
        unit = {"mass": "Jupiter mass"}
        return [{"column_name": c, "datatype": "", "unit": unit.get(c, ""), "description": ""}
                for c in COLS.get(table, [])]

    def simbad_dr3(self, names):
        return {"HOST A": "999"}, "fake"

    def dr3_rows(self, ids):
        out = {}
        for sid in ids:
            ra, dec = self.pos[sid]
            pm = ("10.0", "-5.0") if sid == IDS["E"] else ("0", "0")
            out[sid] = {"source_id": sid, "ra": ra, "dec": dec, "pmra": pm[0], "pmdec": pm[1],
                        "phot_g_mean_mag": "9.0", "ruwe": "1.0", "ref_epoch": "2016.0"}
        return out, "fake"

    def dr4_cone(self, ra, dec, radius_arcsec, schema="gaiadr4"):
        out = []
        for sid, (r0, d0) in self.pos.items():
            if sid == IDS["G"]:
                continue                       # G has no DR4 counterpart
            if abs(r0 - ra) < 1e-6 and abs(d0 - dec) < 1e-6:
                ruwe = "1.4" if sid == IDS["E"] else "1.0"
                pm = ("10.2", "-5.01") if sid == IDS["E"] else ("0", "0")
                out.append({"source_id": sid, "ra": str(r0), "dec": str(d0), "phot_g_mean_mag": "9.02",
                            "ruwe": ruwe, "pmra": pm[0], "pmdec": pm[1], "ref_epoch": "2016.0"})
                out.append({"source_id": "77" + sid[2:], "ra": str(r0), "dec": str(d0 + 1.0 / 3600),
                            "phot_g_mean_mag": "16.0", "ruwe": "1.0", "ref_epoch": "2016.0"})
        return out

    def rows_by_ids(self, table, id_col, ids):
        ids = set(ids)
        return [r for r in self.tables.get(table, []) if r.get(id_col) in ids]

    def count(self, table, where=""):
        return 6000 if where else 7000


def _targets():
    def t(i, group, host, spec="", p=""):
        ra, dec = 0.0, 0.0
        return {"id": i, "group": group, "target": "host " + host, "simbad": "HOST " + host,
                "dr3_source_id": IDS[host], "ra": ra, "dec": dec, "period_spec": spec, "p": p,
                "forecast": "f"}
    return [t("P1", "detectable", "A", "3054.9"), t("P8", "not_detectable", "B", "2540"),
            t("P9", "not_detectable", "C", "418.0"), t("B1", "detectable", "D", "1823.5"),
            t("W1", "w1", "E", "", "0.435"), t("W2", "w2", "E", "", "0.948"),
            t("B3", "b3", "F", "256"), t("P10", "not_detectable", "G", "578.2"),
            {"id": "L1", "group": "l1", "target": "all sky", "simbad": "", "dr3_source_id": "", "ra": "",
             "dec": "", "period_spec": "", "p": "", "forecast": "7500"}]


class EndToEnd(unittest.TestCase):
    def test_not_ready(self):
        with tempfile.TemporaryDirectory() as d:
            res = job.run(FakeIO(ready=False), _targets(), d, log=lambda m: None)
            self.assertEqual(res["stage"], "archive_not_ready")
            self.assertTrue(os.path.exists(os.path.join(d, "results.json")))

    def test_full_run(self):
        with tempfile.TemporaryDirectory() as d:
            res = job.run(FakeIO(), _targets(), d, log=lambda m: None)
            self.assertEqual(res["stage"], "done")
            g = {x["id"]: x for x in res["grades"]}
            self.assertEqual(g["P1"]["status"], "ORBIT_MATCH")
            self.assertEqual(g["P1"]["label"], "RIGHT")
            self.assertEqual(g["P8"]["status"], "NOTHING")
            self.assertEqual(g["P8"]["label"], "RIGHT")
            self.assertEqual(g["P9"]["label"], "WRONG")
            self.assertEqual(g["B1"]["status"], "ACCEL")
            self.assertEqual(g["B1"]["label"], "PARTLY")
            self.assertEqual(g["P10"]["label"], "UNGRADED_NO_COUNTERPART")
            self.assertEqual(g["W1"]["y"], 1)
            self.assertAlmostEqual(g["W1"]["brier"], 0.319225)
            self.assertEqual(g["W2"]["y"], 1)
            self.assertAlmostEqual(g["W2"]["brier"], 0.002704)
            self.assertEqual(g["B3"]["label"], "BH17_SUPPORTED")
            self.assertEqual(g["L1"]["y"], 6000)
            self.assertEqual(g["L1"]["label"], "RIGHT")
            s = res["summary"]
            self.assertEqual(s["detectable"]["graded"], 2)
            self.assertEqual(s["detectable"]["right"], 1)
            self.assertEqual(s["not_detectable"]["graded"], 2)
            self.assertEqual(s["not_detectable"]["right"], 1)
            self.assertEqual(s["binary_all"]["hit_rate"], 0.5)
            self.assertTrue(any("SIMBAD now gives" in p for p in res["problems"]))
            for f in ("results.json", "grades.csv", "solutions.csv", "crossmatch.csv",
                      "schema_tables.csv", "schema_columns.csv", "queries.json"):
                self.assertTrue(os.path.exists(os.path.join(d, f)), f)
            with open(os.path.join(d, "results.json")) as f:
                self.assertEqual(json.load(f)["stage"], "done")
            xm = {x["dr3_source_id"]: x for x in res["crossmatch"]}
            self.assertEqual(xm[IDS["A"]]["dr4_source_id"], IDS["A"])
            self.assertTrue(xm[IDS["A"]]["official"].startswith("dr3_neighbourhood:"))
            self.assertEqual(xm[IDS["G"]]["dr4_source_id"], "")


if __name__ == "__main__":
    unittest.main()

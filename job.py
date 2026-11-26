# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
# Made with AI agents working in MV2's name.
"""Gaia DR4 forecast scorecard: grade the 24 frozen forecasts of PREREGISTRATION.md on release day.

    python job.py --out-dir out [--targets targets.csv] [--schema gaiadr4]

It survives unknown DR4 table names: it dumps the schema, grades what it can, and says what it could not."""
import argparse
import csv
import json
import os
import sys
import time

import scorecard as sc

HERE = os.path.dirname(os.path.abspath(__file__))
NSS_WORDS = ("nss", "orbit", "accel", "exoplanet", "planet", "two_body", "non_single", "substellar",
             "companion", "binary")
SKIP_WORDS = ("epoch", "photometry", "spectra", "rvs_", "xp_", "mean_spectrum", "sampled")
XM_WORDS = ("neighbour", "neighbor", "xmatch", "crossmatch", "best_match")
DEFAULT_DR4_EPOCH = 2017.5
MJUP_IN_MSUN = 0.000954588
MJUP_IN_MEARTH = 317.83


def load_targets(path):
    with open(path, encoding="utf-8") as f:
        return sc.read_csv_text(f.read())


def short(t):
    return t.split(".", 1)[-1].lower()


def write_csv(path, rows, cols=None):
    cols = cols or sorted({k for r in rows for k in r})
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in cols})


def nss_like(table):
    t = short(table)
    return any(w in t for w in NSS_WORDS) and not any(w in t for w in SKIP_WORDS)


def xm_like(table):
    t = short(table)
    return "dr3" in t and any(w in t for w in XM_WORDS)


def mass_threshold(unit):
    u = (unit or "").lower()
    if "jup" in u:
        return 13.0
    if "sun" in u or "sol" in u:
        return 13.0 * MJUP_IN_MSUN
    if "earth" in u or "geo" in u:
        return 13.0 * MJUP_IN_MEARTH
    return None


def run(io, targets, out_dir, schema="gaiadr4", log=print):
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()
    res = {"version": sc.VERSION, "preregistration": "PREREGISTRATION.md (frozen 2026-11-25)",
           "schema": schema, "problems": []}

    def problem(msg):
        log("PROBLEM " + msg)
        res["problems"].append(msg)

    ok, msg = io.preflight(schema)
    res["preflight"] = msg
    if not ok:
        res["stage"] = "archive_not_ready"
        log("archive not ready: %s" % msg)
        with open(os.path.join(out_dir, "results.json"), "w") as f:
            json.dump(res, f, indent=1)
        return res

    # 1. schema dump
    try:
        tables = io.list_tables(schema)
    except Exception as e:
        tables = []
        problem("list_tables failed: %s" % e)
    write_csv(os.path.join(out_dir, "schema_tables.csv"), tables, ["schema_name", "table_name", "description"])
    names = [t.get("table_name", "") for t in tables]
    names = [n if "." in n else "%s.%s" % (schema, n) for n in names if n]
    nss_tables = [t for t in names if nss_like(t)]
    xm_tables = [t for t in names if xm_like(t)]
    log("DR4 tables: %d; NSS-like: %s; crossmatch-like: %s" % (len(names), nss_tables, xm_tables))
    columns, col_rows = {}, []
    for t in nss_tables + xm_tables:
        cr = []
        for name in (t, short(t)):
            try:
                cr = io.list_columns(name)
            except Exception as e:
                problem("list_columns %s failed: %s" % (name, e))
                cr = []
            if cr:
                break
        columns[t] = {c["column_name"].lower(): c for c in cr if c.get("column_name")}
        for c in cr:
            col_rows.append(dict(c, table_name=t))
    write_csv(os.path.join(out_dir, "schema_columns.csv"), col_rows,
              ["table_name", "column_name", "datatype", "unit", "description"])
    res["nss_tables"], res["crossmatch_tables"] = nss_tables, xm_tables

    # 2. identifiers: frozen, checked against SIMBAD
    hosts = {}
    for t in targets:
        sid = sc.id_str(t.get("dr3_source_id"))
        if sid:
            hosts.setdefault(sid, t)
    try:
        sim, src = io.simbad_dr3([t.get("simbad") for t in hosts.values()])
        res["simbad"] = src
        for sid, t in hosts.items():
            now = sim.get(t.get("simbad"))
            if now and now != sid:
                problem("SIMBAD now gives Gaia DR3 %s for %s (frozen %s); frozen id kept" % (now, t["simbad"], sid))
    except Exception as e:
        res["simbad"] = "failed: %s" % e

    # 3. DR3 rows
    try:
        dr3, dr3_src = io.dr3_rows(list(hosts))
    except Exception as e:
        dr3, dr3_src = {}, "failed: %s" % e
        problem("DR3 rows failed: %s" % e)
    res["dr3_source"] = dr3_src

    # 4. DR4 counterparts
    xm, dr4_of, dr4_row = [], {}, {}
    dr4_epoch = None
    for sid, t in hosts.items():
        d3 = dr3.get(sid)
        if d3 is None:
            d3 = {"source_id": sid, "ra": t["ra"], "dec": t["dec"], "pmra": None, "pmdec": None,
                  "phot_g_mean_mag": None}
            problem("no DR3 row for %s; SIMBAD position used without proper motion" % t["simbad"])
        try:
            cands = io.dr4_cone(float(d3["ra"]), float(d3["dec"]), 5.0, schema)
        except Exception as e:
            problem("DR4 cone failed for %s: %s" % (t["simbad"], e))
            cands = []
        if dr4_epoch is None:
            for c in cands:
                if sc.to_float(c.get("ref_epoch")):
                    dr4_epoch = sc.to_float(c.get("ref_epoch"))
                    break
        ep = dr4_epoch or DEFAULT_DR4_EPOCH
        best, info = sc.crossmatch_pick(d3, cands, ep, sc.to_float(d3.get("ref_epoch")) or 2016.0)
        row = {"simbad": t["simbad"], "dr3_source_id": sid, "n_cone": len(cands),
               "dr4_source_id": sc.id_str(best.get("source_id")) if best else "",
               "sep_arcsec": info["sep_arcsec"], "n_candidates_ok": info["n_candidates_ok"],
               "same_source_id": info["same_source_id"], "dr4_epoch_used": ep, "official": ""}
        # official crossmatch, if the archive has one: it wins, disagreements are reported
        for xt in xm_tables:
            cols = columns.get(xt, {})
            c3 = next((c for c in cols if "dr3" in c and "source_id" in c), None)
            c4 = next((c for c in cols if "dr4" in c and "source_id" in c), None) or (
                "source_id" if "source_id" in cols and c3 != "source_id" else None)
            if not (c3 and c4):
                continue
            try:
                xr = io.rows_by_ids(xt, c3, [sid])
            except Exception as e:
                problem("crossmatch table %s failed: %s" % (xt, e))
                continue
            if xr:
                off = sc.id_str(xr[0].get(c4))
                row["official"] = "%s:%s" % (short(xt), off)
                if off and off != row["dr4_source_id"]:
                    problem("official crossmatch %s gives %s for %s, position rule gave %s; official used"
                            % (xt, off, t["simbad"], row["dr4_source_id"]))
                    row["dr4_source_id"] = off
                    best = next((c for c in cands if sc.id_str(c.get("source_id")) == off), None)
                break
        xm.append(row)
        if row["dr4_source_id"]:
            dr4_of[sid] = row["dr4_source_id"]
            if best is None:
                try:
                    got = io.rows_by_ids("%s.gaia_source" % schema, "source_id", [row["dr4_source_id"]])
                    best = got[0] if got else None
                except Exception as e:
                    problem("DR4 row fetch failed for %s: %s" % (t["simbad"], e))
            dr4_row[sid] = best or {}
    write_csv(os.path.join(out_dir, "crossmatch.csv"), xm)

    # 5. non-single-star rows for the counterparts
    sols = {sid: [] for sid in hosts}
    sol_rows = []
    by4 = {v: k for k, v in dr4_of.items()}
    for nt in nss_tables:
        if "source_id" not in columns.get(nt, {}):
            continue
        try:
            rows = io.rows_by_ids(nt, "source_id", list(by4))
        except Exception as e:
            problem("NSS query %s failed: %s" % (nt, e))
            continue
        for r in rows:
            sid = by4.get(sc.id_str(r.get("source_id")))
            if sid is None:
                continue
            s = sc.solution_from_row(r, nt)
            sols[sid].append(s)
            sol_rows.append(dict(s, simbad=hosts[sid]["simbad"], dr4_source_id=sc.id_str(r.get("source_id"))))
    write_csv(os.path.join(out_dir, "solutions.csv"), sol_rows,
              ["simbad", "dr4_source_id", "table", "solution_type", "kind", "period", "inclination"])

    # 6. grades
    grades = []
    for t in targets:
        g = {"id": t["id"], "group": t["group"], "target": t["target"], "forecast": t["forecast"],
             "period_spec": t.get("period_spec", ""), "status": "", "label": "", "y": "", "brier": "",
             "ruwe_dr4": "", "detail": ""}
        sid = sc.id_str(t.get("dr3_source_id"))
        have = sid in dr4_of
        if have:
            g["ruwe_dr4"] = dr4_row.get(sid, {}).get("ruwe", "")
        if t["group"] in (sc.DETECTABLE, sc.NOT_DETECTABLE):
            status = sc.host_status(sols.get(sid, []), t["period_spec"]) if have else None
            g["status"] = status or ""
            g["label"] = sc.label_for(t["group"], status)
            orb = [s for s in sols.get(sid, []) if s["kind"] != "other"]
            g["detail"] = "; ".join("%s %s P=%s" % (s["kind"], s["solution_type"], s["period"]) for s in orb)
            if g["label"] == "WRONG":
                g["detail"] += " | WRONG means DR4 did not show it, not that the planet is absent"
        elif t["group"] == "w1":
            y = sc.w1_outcome(g["ruwe_dr4"]) if have else None
            g["y"] = "" if y is None else y
            g["brier"] = "" if y is None else round(sc.brier(t["p"], y), 6)
            g["label"] = "GRADED" if y is not None else "UNGRADED"
        elif t["group"] == "w2":
            d3, d4 = dr3.get(sid, {}), dr4_row.get(sid, {})
            y, dra, ddec = sc.w2_outcome(d4.get("pmra"), d4.get("pmdec"), d3.get("pmra"), d3.get("pmdec"))
            g["y"] = "" if y is None else y
            g["brier"] = "" if y is None else round(sc.brier(t["p"], y), 6)
            g["label"] = "GRADED" if y is not None else "UNGRADED"
            g["detail"] = "" if y is None else "dmu_ra=%.4f dmu_dec=%.4f mas/yr" % (dra, ddec)
        elif t["group"] == "b3":
            status = sc.host_status(sols.get(sid, []), t["period_spec"]) if have else None
            inc = None
            for s in sols.get(sid, []):
                if s["kind"] == "orbit" and sc.period_match(s["period"], t["period_spec"]):
                    inc = s["inclination"]
                    break
            g["status"] = status or ""
            g["label"] = sc.b3_label(status, inc)
            g["detail"] = "" if inc is None else "inclination %.2f deg" % inc
        elif t["group"] == "l1":
            n, how = None, "no planet-like DR4 table found"
            for nt in nss_tables:
                if not any(w in short(nt) for w in ("exoplanet", "planet", "substellar")):
                    continue
                try:
                    total = io.count(nt)
                except Exception as e:
                    problem("count %s failed: %s" % (nt, e))
                    continue
                n, how = total, "%s rows: %s" % (short(nt), total)
                cols = columns.get(nt, {})
                mc = next((c for c in cols if "mass" in c and "error" not in c and "err" not in c), None)
                thr = mass_threshold(cols.get(mc, {}).get("unit")) if mc else None
                if mc and thr:
                    try:
                        nm = io.count(nt, "%s < %.6g" % (mc, thr))
                        n, how = nm, how + "; %s < 13 MJup: %s" % (mc, nm)
                    except Exception as e:
                        problem("mass count %s failed: %s" % (nt, e))
                break
            g["y"] = "" if n is None else n
            g["label"] = sc.l1_label(n)
            g["detail"] = how + " | PROVISIONAL: the frozen rule uses ESA's stated number if the release documentation gives one"
        grades.append(g)
    write_csv(os.path.join(out_dir, "grades.csv"), grades,
              ["id", "group", "target", "forecast", "period_spec", "status", "label", "y", "brier", "ruwe_dr4", "detail"])

    res["summary"] = sc.summarize(grades)
    res["grades"] = grades
    res["crossmatch"] = xm
    res["stage"] = "done"
    res["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(out_dir, "results.json"), "w") as f:
        json.dump(res, f, indent=1, default=str)
    queries = getattr(io, "QUERY_LOG", [])
    with open(os.path.join(out_dir, "queries.json"), "w") as f:
        json.dump(queries, f, indent=1)
    log("OUTCOME %s" % json.dumps(res["summary"]))
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--targets", default=os.path.join(HERE, "targets.csv"))
    ap.add_argument("--schema", default="gaiadr4")
    a = ap.parse_args(argv)
    import io_dr4score
    res = run(io_dr4score, load_targets(a.targets), a.out_dir, a.schema,
              log=lambda m: print("[%7.1fs] %s" % (time.time() - T0, m), flush=True))
    return 0 if res.get("stage") == "done" else 3


T0 = time.time()
if __name__ == "__main__":
    sys.exit(main())

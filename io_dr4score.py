# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
# Made with AI agents working in MV2's name.
"""Network I/O for the Gaia DR4 forecast scorecard (standard library only, CSV answers).

Gaia archive TAP for DR4 (and DR3 when it is served again), SIMBAD TAP (two mirrors) for the frozen
host identifiers, VizieR TAP (I/355/gaiadr3) as the DR3 fallback. No analysis here; rules are in
scorecard.py. Every query string is logged by job.py."""
import time
import urllib.parse
import urllib.request

import scorecard as sc

VERSION = "0.1.0"
UA = {"User-Agent": "gaia-dr4-forecast-scorecard/0.1 (open research code; Apache-2.0)"}
GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
SIMBAD_TAPS = ("https://simbad.cds.unistra.fr/simbad/sim-tap/sync",
               "https://simbad.cfa.harvard.edu/simbad/sim-tap/sync")
VIZIER_TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"

QUERY_LOG = []


def fetch(url, timeout=300):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def tap(base, query, tries=3, fmt_key="FORMAT"):
    params = {"REQUEST": "doQuery", "LANG": "ADQL", fmt_key: "csv", "QUERY": query}
    url = base + "?" + urllib.parse.urlencode(params)
    QUERY_LOG.append({"service": base, "query": query})
    last = None
    for k in range(tries):
        try:
            text = fetch(url).decode("utf-8", "replace")
            if text.lstrip().lower().startswith("<"):
                raise RuntimeError("non-CSV answer: %s" % text[:300])
            return sc.read_csv_text(text)
        except Exception as e:
            last = e
            time.sleep(5 * (k + 1))
    raise RuntimeError("TAP failed (%s): %s" % (base, str(last)[:300]))


def gaia(query, tries=3):
    return tap(GAIA_TAP, query, tries)


# ---------------------------------------------------------------- Gaia archive
def preflight(schema="gaiadr4"):
    try:
        rows = gaia("SELECT TOP 1 source_id FROM %s.gaia_source" % schema, tries=2)
        return bool(rows) and "source_id" in rows[0], "%d row(s)" % len(rows)
    except Exception as e:
        return False, str(e)[:300]


def list_tables(schema="gaiadr4"):
    # descriptions are left out: free text with line breaks would upset the line-based CSV reader
    return gaia("SELECT schema_name, table_name FROM TAP_SCHEMA.tables "
                "WHERE schema_name = '%s'" % schema)


def list_columns(table):
    return gaia("SELECT column_name, datatype, unit FROM TAP_SCHEMA.columns "
                "WHERE table_name = '%s'" % table)


def dr4_cone(ra, dec, radius_arcsec, schema="gaiadr4"):
    q = ("SELECT * FROM %s.gaia_source WHERE 1=CONTAINS(POINT('ICRS', ra, dec), "
         "CIRCLE('ICRS', %.8f, %.8f, %.8f))" % (schema, ra, dec, radius_arcsec / 3600.0))
    return gaia(q)


def rows_by_ids(table, id_col, ids):
    ids = [sc.id_str(i) for i in ids if sc.id_str(i)]
    if not ids:
        return []
    return gaia("SELECT * FROM %s WHERE %s IN (%s)" % (table, id_col, ", ".join(ids)))


def count(table, where=""):
    rows = gaia("SELECT COUNT(*) AS n FROM %s%s" % (table, (" WHERE " + where) if where else ""))
    return int(float(rows[0]["n"])) if rows else None


# ---------------------------------------------------------------- DR3 (archive, then VizieR)
DR3_COLS = "source_id, ra, dec, pmra, pmdec, phot_g_mean_mag, ruwe, ref_epoch"


def dr3_rows(ids):
    """{source_id: row} with ra, dec, pmra, pmdec, phot_g_mean_mag, ruwe (epoch 2016.0)."""
    ids = [sc.id_str(i) for i in ids if sc.id_str(i)]
    if not ids:
        return {}, "no ids"
    try:
        rows = gaia("SELECT %s FROM gaiadr3.gaia_source WHERE source_id IN (%s)"
                    % (DR3_COLS, ", ".join(ids)))
        if rows:
            return {sc.id_str(r["source_id"]): r for r in rows}, "gaia archive gaiadr3.gaia_source"
    except Exception:
        pass
    rows = tap(VIZIER_TAP, 'SELECT Source, RA_ICRS, DE_ICRS, pmRA, pmDE, Gmag, RUWE FROM "I/355/gaiadr3" '
               "WHERE Source IN (%s)" % ", ".join(ids), fmt_key="RESPONSEFORMAT")
    out = {}
    for r in rows:
        sid = sc.id_str(r.get("source"))
        out[sid] = {"source_id": sid, "ra": r.get("ra_icrs"), "dec": r.get("de_icrs"),
                    "pmra": r.get("pmra"), "pmdec": r.get("pmde"), "phot_g_mean_mag": r.get("gmag"),
                    "ruwe": r.get("ruwe"), "ref_epoch": "2016.0"}
    return out, "VizieR I/355/gaiadr3"


# ---------------------------------------------------------------- SIMBAD
def simbad_dr3(names):
    """{simbad name: 'Gaia DR3 <id>' digits} from the first SIMBAD mirror that answers."""
    names = [n for n in names if n]
    if not names:
        return {}, "no names"
    inlist = ", ".join("'%s'" % n.replace("'", "''") for n in names)
    q = ("SELECT id1.id AS name, id2.id AS gaia FROM ident AS id1 JOIN ident AS id2 "
         "ON id1.oidref = id2.oidref WHERE id1.id IN (%s) AND id2.id LIKE 'Gaia DR3%%'" % inlist)
    last = None
    for base in SIMBAD_TAPS:
        try:
            rows = tap(base, q, tries=2)
            return {" ".join(r["name"].split()): sc.id_str(r["gaia"]) for r in rows}, base
        except Exception as e:
            last = e
    return {}, "SIMBAD failed: %s" % str(last)[:200]

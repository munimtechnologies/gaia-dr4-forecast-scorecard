# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
# Made with AI agents working in MV2's name.
"""Grading rules of the Gaia DR4 forecast scorecard (PREREGISTRATION.md section 4). Pure functions,
standard library only; all network I/O lives in io_dr4score.py."""
import csv
import io
import math

VERSION = "0.1.0"

RUWE_CUT = 1.15            # W1
W1_P = 0.435
W2_P = 0.948
W2_SIG_RA = 0.030          # mas/yr, the paper's combined errors
W2_SIG_DEC = 0.025
W2_NSIG = 3.0              # frozen choice
PERIOD_TOL = 0.20          # ORBIT_MATCH: within 20% of the forecast period
XM_RADIUS_ARCSEC = 2.0
XM_DMAG = 0.5
BH17_INC = 10.9
BH17_TOL = 5.0
NOT_BH17_MIN = 30.0
L1_LO, L1_HI = 3300, 11700

DETECTABLE = "detectable"
NOT_DETECTABLE = "not_detectable"


# ---------------------------------------------------------------- small helpers
def to_float(v):
    if v is None:
        return None
    if isinstance(v, (int, float)):
        return None if (isinstance(v, float) and math.isnan(v)) else float(v)
    s = str(v).strip()
    if s == "" or s.lower() in ("null", "nan", "none", "--"):
        return None
    try:
        x = float(s)
    except ValueError:
        return None
    return None if math.isnan(x) else x


def id_str(v):
    """Gaia source_id as a digit string (never via float: 19-digit ids lose precision)."""
    s = str(v if v is not None else "").strip()
    if s.lower().startswith("gaia dr3 ") or s.lower().startswith("gaia dr4 "):
        s = s[9:].strip()
    return s if s.isdigit() else ""


def read_csv_text(text):
    """Rows of a CSV text as dicts with lower-case keys; comment lines (#) skipped."""
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        return []
    rd = csv.DictReader(io.StringIO("\n".join(lines)))
    out = []
    for r in rd:
        out.append({(k or "").strip().lower(): (v.strip() if isinstance(v, str) else v)
                    for k, v in r.items()})
    return out


def parse_periods(spec):
    """'3054.9' -> [(3054.9, 3054.9)]; '1700;3000' -> two; '200-300' -> [(200, 300)]; '' -> []."""
    out = []
    for part in str(spec or "").split(";"):
        part = part.strip()
        if not part:
            continue
        if "-" in part[1:]:
            i = part.index("-", 1)
            lo, hi = float(part[:i]), float(part[i + 1:])
            out.append((min(lo, hi), max(lo, hi)))
        else:
            x = float(part)
            out.append((x, x))
    return out


def period_match(p, spec, tol=PERIOD_TOL):
    """True if period p (days) is within tol of any forecast period (ranges widened by tol at both ends)."""
    p = to_float(p)
    if p is None:
        return False
    return any(lo * (1 - tol) <= p <= hi * (1 + tol) for lo, hi in parse_periods(spec))


# ---------------------------------------------------------------- solutions
def classify_solution(solution_type, table_name="", has_period=False):
    """'orbit' (an orbital solution that includes astrometry), 'accel' (astrometric acceleration or
    similar non-linear motion) or 'other' (spectroscopic-only, eclipsing, trend, ...)."""
    st = (solution_type or "").strip().lower()
    tn = (table_name or "").strip().lower()
    if st:
        if "orbital" in st or "astrospectro" in st or "astrometricorbit" in st or "astro_orbit" in st:
            return "orbit"
        if "acceleration" in st:
            return "accel"
    if any(w in tn for w in ("exoplanet", "planet", "substellar")) and has_period:
        return "orbit"
    if not st and "accel" in tn and "spectro" not in tn:
        return "accel"
    return "other"


def host_status(solutions, period_spec):
    """solutions: list of {'kind': orbit|accel|other, 'period': days or None}.
    ORBIT_MATCH > ORBIT > ACCEL > NOTHING (PREREGISTRATION section 4)."""
    orbits = [s for s in solutions if s.get("kind") == "orbit"]
    if any(period_match(s.get("period"), period_spec) for s in orbits):
        return "ORBIT_MATCH"
    if orbits:
        return "ORBIT"
    if any(s.get("kind") == "accel" for s in solutions):
        return "ACCEL"
    return "NOTHING"


def label_detectable(status):
    return {"ORBIT_MATCH": "RIGHT", "ORBIT": "PARTLY", "ACCEL": "PARTLY", "NOTHING": "WRONG"}[status]


def label_not_detectable(status):
    return "WRONG" if status == "ORBIT_MATCH" else "RIGHT"


def label_for(group, status):
    if status is None:
        return "UNGRADED_NO_COUNTERPART"
    if group == DETECTABLE:
        return label_detectable(status)
    if group == NOT_DETECTABLE:
        return label_not_detectable(status)
    raise ValueError(group)


# ---------------------------------------------------------------- probabilistic forecasts
def brier(p, y):
    return (float(p) - float(y)) ** 2


def w1_outcome(ruwe):
    r = to_float(ruwe)
    if r is None:
        return None
    return 1 if r > RUWE_CUT else 0


def w2_outcome(pmra4, pmdec4, pmra3, pmdec3):
    """(y, dmu_ra, dmu_dec); y = 1 if either axis changed by >= 3 sigma of the paper's errors."""
    v = [to_float(x) for x in (pmra4, pmdec4, pmra3, pmdec3)]
    if any(x is None for x in v):
        return None, None, None
    dra, ddec = v[0] - v[2], v[1] - v[3]
    y = 1 if (abs(dra) / W2_SIG_RA >= W2_NSIG or abs(ddec) / W2_SIG_DEC >= W2_NSIG) else 0
    return y, dra, ddec


# ---------------------------------------------------------------- B3 inclination
def thiele_innes_forward(a, inc_deg, omega_deg, Omega_deg):
    """(A, B, F, G) from Campbell elements; used by tests and as documentation."""
    i, w, W = (math.radians(x) for x in (inc_deg, omega_deg, Omega_deg))
    A = a * (math.cos(w) * math.cos(W) - math.sin(w) * math.sin(W) * math.cos(i))
    B = a * (math.cos(w) * math.sin(W) + math.sin(w) * math.cos(W) * math.cos(i))
    F = a * (-math.sin(w) * math.cos(W) - math.cos(w) * math.sin(W) * math.cos(i))
    G = a * (-math.sin(w) * math.sin(W) + math.cos(w) * math.cos(W) * math.cos(i))
    return A, B, F, G


def inclination_from_thiele_innes(A, B, F, G):
    """Inclination (deg, 0-180) from the Thiele-Innes constants: with p = a(1 + cos i) and
    m = a(1 - cos i), cos i = (p - m) / (p + m)."""
    v = [to_float(x) for x in (A, B, F, G)]
    if any(x is None for x in v):
        return None
    A, B, F, G = v
    p = math.hypot(A + G, B - F)
    m = math.hypot(A - G, B + F)
    if p + m == 0:
        return None
    c = max(-1.0, min(1.0, (p - m) / (p + m)))
    return math.degrees(math.acos(c))


def b3_label(status, inc_deg):
    if status != "ORBIT_MATCH":
        return "NO_ORBIT" if status in ("NOTHING", "ACCEL", None) else "NO_MATCHING_ORBIT"
    i = to_float(inc_deg)
    if i is None:
        return "MATCH_NO_INCLINATION"
    j = min(i, 180.0 - i)
    if abs(j - BH17_INC) <= BH17_TOL:
        return "BH17_SUPPORTED"
    if j >= NOT_BH17_MIN:
        return "NOT_BH17"
    return "INTERMEDIATE"


def l1_label(n):
    n = to_float(n)
    if n is None:
        return "UNGRADABLE"
    if n < L1_LO:
        return "LOW"
    if n > L1_HI:
        return "HIGH"
    return "RIGHT"


# ---------------------------------------------------------------- crossmatch
def sep_arcsec(ra1, dec1, ra2, dec2):
    r = math.radians
    a = (math.sin(r(dec2 - dec1) / 2) ** 2
         + math.cos(r(dec1)) * math.cos(r(dec2)) * math.sin(r(ra2 - ra1) / 2) ** 2)
    return math.degrees(2 * math.asin(min(1.0, math.sqrt(a)))) * 3600.0


def propagate(ra, dec, pmra, pmdec, dt_yr):
    """Linear propagation; pmra includes cos(dec) (Gaia convention), mas/yr."""
    pmra = to_float(pmra) or 0.0
    pmdec = to_float(pmdec) or 0.0
    cosd = max(1e-9, math.cos(math.radians(dec)))
    return (ra + pmra * dt_yr / 3.6e6 / cosd, dec + pmdec * dt_yr / 3.6e6)


def crossmatch_pick(dr3, candidates, dr4_epoch, dr3_epoch=2016.0):
    """Nearest DR4 candidate within 2 arcsec of the propagated DR3 position with |dG| < 0.5.
    dr3: dict with ra, dec, pmra, pmdec, phot_g_mean_mag. Returns (row or None, info dict)."""
    ra3, dec3 = to_float(dr3.get("ra")), to_float(dr3.get("dec"))
    g3 = to_float(dr3.get("phot_g_mean_mag"))
    ra_p, dec_p = propagate(ra3, dec3, dr3.get("pmra"), dr3.get("pmdec"), dr4_epoch - dr3_epoch)
    best, best_sep, n_close = None, None, 0
    for c in candidates:
        ra4, dec4 = to_float(c.get("ra")), to_float(c.get("dec"))
        if ra4 is None or dec4 is None:
            continue
        s = sep_arcsec(ra_p, dec_p, ra4, dec4)
        if s > XM_RADIUS_ARCSEC:
            continue
        g4 = to_float(c.get("phot_g_mean_mag"))
        if g3 is not None and g4 is not None and abs(g4 - g3) >= XM_DMAG:
            continue
        n_close += 1
        if best_sep is None or s < best_sep:
            best, best_sep = c, s
    info = {"sep_arcsec": best_sep, "n_candidates_ok": n_close,
            "same_source_id": (best is not None
                               and id_str(best.get("source_id")) == id_str(dr3.get("source_id")))}
    return best, info


# ---------------------------------------------------------------- column finders
def find_col(row_or_cols, *names, contains=None):
    """First column (lower-case) equal to one of names, else the first containing `contains`."""
    cols = list(row_or_cols.keys()) if isinstance(row_or_cols, dict) else list(row_or_cols)
    low = [c.lower() for c in cols]
    for n in names:
        if n in low:
            return n
    if contains:
        for c in low:
            if contains in c:
                return c
    return None


def solution_from_row(row, table_name):
    """{'kind', 'period', 'solution_type', 'inclination', 'table'} from one NSS-like row."""
    st_col = find_col(row, "nss_solution_type", "solution_type", contains="solution_type")
    p_col = find_col(row, "period", "period_days", "orbital_period")
    st = row.get(st_col) if st_col else None
    per = to_float(row.get(p_col)) if p_col else None
    inc = None
    i_col = find_col(row, "inclination", "inc")
    if i_col:
        inc = to_float(row.get(i_col))
    if inc is None:
        inc = inclination_from_thiele_innes(row.get("a_thiele_innes"), row.get("b_thiele_innes"),
                                            row.get("f_thiele_innes"), row.get("g_thiele_innes"))
    return {"table": table_name, "solution_type": st or "", "period": per, "inclination": inc,
            "kind": classify_solution(st, table_name, per is not None)}


# ---------------------------------------------------------------- summary
def summarize(grades):
    """grades: list of dicts with 'group' and 'label'. Hit rates over graded binary forecasts."""
    out = {}
    for g in (DETECTABLE, NOT_DETECTABLE):
        rows = [x for x in grades if x["group"] == g and x["label"] in ("RIGHT", "PARTLY", "WRONG")]
        n_right = sum(1 for x in rows if x["label"] == "RIGHT")
        out[g] = {"graded": len(rows), "right": n_right,
                  "partly": sum(1 for x in rows if x["label"] == "PARTLY"),
                  "wrong": sum(1 for x in rows if x["label"] == "WRONG"),
                  "hit_rate": (n_right / len(rows)) if rows else None}
    allg = out[DETECTABLE]["graded"] + out[NOT_DETECTABLE]["graded"]
    allr = out[DETECTABLE]["right"] + out[NOT_DETECTABLE]["right"]
    out["binary_all"] = {"graded": allg, "right": allr, "hit_rate": (allr / allg) if allg else None}
    return out

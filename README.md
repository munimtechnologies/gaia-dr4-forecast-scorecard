# gaia-dr4-forecast-scorecard

A pre-registered, day-one scorecard of published forecasts for Gaia DR4 (expected 2026-12-02): which
named stars would show their companions, and how many giant planets DR4 would find.

- **Status:** pre-registered 2026-11-25, before any DR4 data exist. Rules: [PREREGISTRATION.md](PREREGISTRATION.md).
  Code v0.1.0 committed 2026-11-26, also before release.
- **Forecasts graded:** 24 forecasts (20 yes/no, 2 probabilistic, 1 inclination test, 1 yield) from four papers
  (Wallace, Wittenmyer & Horner, arXiv:2610.02783; Parc et al., arXiv:2610.25152; Baycroft, Triaud & Sahlmann,
  arXiv:2603.02198; Lammers & Winn, AJ 171, 18).

## Code (standard library only)

- `scorecard.py`: the frozen rules of PREREGISTRATION.md section 4 (ORBIT / ORBIT_MATCH / ACCEL / NOTHING,
  RIGHT / PARTLY / WRONG, Brier scores, the HD 202206 inclination test with inclinations from Thiele-Innes
  constants, the yield label, the DR3-to-DR4 crossmatch).
- `io_dr4score.py`: Gaia archive TAP (DR4; DR3 when served), SIMBAD TAP, VizieR TAP (I/355/gaiadr3) fallback.
- `job.py`: `python job.py --out-dir out`. It dumps the DR4 schema first (`schema_tables.csv`,
  `schema_columns.csv`) and uses every non-single-star-like table it finds, so unknown DR4 table names do not stop
  it. Outputs: `grades.csv`, `solutions.csv`, `crossmatch.csv`, `results.json`, `queries.json`.
- `targets.csv`: the 24 forecasts as frozen, with Gaia DR3 identifiers from SIMBAD (re-checked at run time).
- `test_scorecard.py`: `python -m unittest -v test_scorecard` (fake archive, no network).
- `dr4day/job.template.json`: the release-day job.

**Readings fixed before release (not changes to the rules):** where a paper gives two candidate periods
(WASP-89: about 1,700 or 3,000 d) or a range (WASP-79: 200-300 d), ORBIT_MATCH means within 20% of any of them.
BEBOP-4 is HD 111605 (Triaud et al. 2025, arXiv:2510.25329; companion period 1,823.5 d around a 72-d binary), so
an orbit at the 72-d binary period is reported, not counted. The L1 count the job prints is provisional: the
frozen rule uses ESA's own stated number of DR4 exoplanet candidates when the release documentation gives one.

By MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name; no human expert has
checked it yet.

Licence: Apache-2.0 (code); text and data CC BY 4.0.

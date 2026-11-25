# Pre-registration: a day-one scorecard of published Gaia DR4 forecasts

Frozen 2026-11-25 (UTC), seven days before Gaia DR4 (expected 2026-12-02). No DR4 data exist yet, and
the Gaia archive is closed until the release. MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents
working in MV2's name; no human expert has checked it yet. Licence: Apache-2.0 (code), CC BY 4.0 (text, data).

## 1. Question

Several 2026 papers made dated, checkable forecasts of what Gaia DR4 would show for named stars and for
the whole planet yield. Once DR4 is out, each forecast is right or wrong, but no single paper will grade
the others. We will grade all forecasts frozen here, with the same rules, as soon as DR4 is public.

The interesting result, if it comes, is a measured accuracy for the forecasting methods everyone is now
using to plan DR4 follow-up (Lammers & Winn's SNR1 criterion, Baycroft et al.'s SNR criterion,
posterior-predictive RUWE forecasts). A forecast method that over- or under-predicts detections by a
large factor would change how the DR4 planet lists should be read.

## 2. Forecasts (frozen; quoted from the papers)

| ID | Source | Target | Forecast for DR4 |
|---|---|---|---|
| W1 | Wallace, Wittenmyer & Horner 2026, arXiv:2610.02783 (PASA) | HD 65216 | P(RUWE > 1.15) = 0.435 (43.5% of posterior samples) |
| W2 | same | HD 65216 | P(detectable DR3-to-DR4 proper-motion change) = 0.948; combined errors 0.030 (RA) and 0.025 (Dec) mas/yr |
| P1-P7 | Parc et al. 2026, arXiv:2610.25152 (A&A), Lammers & Winn SNR1 method | WASP-8 c, WASP-34 c, WASP-53 c, WASP-68 c, WASP-81 c, WASP-99 c, WASP-132 d | detectable in DR4 (SNR1 2.2-11.5; WASP-34 c 2.23 one-companion model) |
| P8-P13 | same | WASP-22, WASP-41 c, WASP-47 c, WASP-67, WASP-134 c, WASP-185 c | not detectable in DR4 |
| P14-P17 | same (ambiguous RV signals) | WASP-89, WASP-7, WASP-79, WASP-84 | not detectable in DR4 (WASP-89 SNR1 = 1.09, detectable only in DR5) |
| P18 | same (trend only) | WASP-171 | "could expect the detection" (weak forecast: detectable) |
| B1 | Baycroft, Triaud & Sahlmann 2026, arXiv:2603.02198 (MNRAS) | BEBOP-4 (circumbinary brown dwarf, about 20 MJup) | SNR about 40, above their threshold of 20: detectable |
| B2 | same | Kepler-16 b | SNR about 1: not detectable |
| B3 | same | HD 202206 (inner/outer companions) | SNR about 80/70 if edge-on; about 400/250 if the Benedict & Harrison (2017) inclinations (10.9 and 7.7 deg) are right |
| L1 | Lammers & Winn 2026, AJ 171, 18 | all sky | 7,500 +- 2,100 giant planets detectable in DR4 |

Our own forecast test of Proxima c (munimtechnologies/proxima-c-dr4) is graded in its own repository
and is not part of this scorecard.

## 3. Data (all public on release day)

- Gaia DR4 through the ESA Gaia archive TAP service (gaia_source; the non-single-star tables, whatever
  their DR4 names; any exoplanet-candidate table). Gaia DR3 gaia_source for W2.
- Host positions and Gaia DR3 identifiers from SIMBAD (TAP). DR4 counterpart: the nearest DR4 source
  within 2 arcsec of the DR3 position propagated to the DR4 epoch, with |G_DR4 - G_DR3| < 0.5 mag; if
  the archive provides an official DR3-DR4 crossmatch, that is used instead and any disagreement reported.
- ESA's own release documentation for the number of exoplanet candidates (L1).

## 4. Grading rules (fixed now)

Definitions for a host in DR4:
- **ORBIT**: the host has an astrometric orbital solution (any orbital solution type that includes
  astrometry) in the DR4 non-single-star or exoplanet-candidate tables.
- **ORBIT_MATCH**: ORBIT, and the solution period is within 20% of the forecast companion's period
  (RV period from the forecasting paper; for B3 the inner companion, about 256 d).
- **ACCEL**: no ORBIT, but an astrometric acceleration (or similar non-linear motion) solution.
- **NOTHING**: none of these. RUWE is reported for every host.

Per forecast:
- "Detectable" forecasts (P1-P7, P18, B1) are **RIGHT** if ORBIT_MATCH, **PARTLY** if ORBIT without
  match or ACCEL, **WRONG** if NOTHING.
- "Not detectable" forecasts (P8-P17, B2) are **RIGHT** if there is no ORBIT_MATCH, **WRONG** if
  ORBIT_MATCH. An orbit at an unrelated period (for example a stellar companion) is reported, not counted.
- W1: outcome y = 1 if DR4 RUWE > 1.15, else 0. Brier score (0.435 - y)^2.
- W2: dmu = mu_DR4 - mu_DR3 per axis; y = 1 if |dmu_ra|/0.030 >= 3 or |dmu_dec|/0.025 >= 3 (mas/yr;
  the 3-sigma level is our choice, frozen here; the paper's errors). Brier score (0.948 - y)^2.
- B3: if HD 202206 has ORBIT_MATCH with an inclination i, then **BH17_SUPPORTED** if min(i, 180 - i)
  is within 5 deg of 10.9; **NOT_BH17** if min(i, 180 - i) >= 30; else **INTERMEDIATE**. Without an
  orbit: **NO_ORBIT** (itself evidence against SNR 400 if the pipeline solved comparable systems).
- L1: N = ESA's stated number of DR4 exoplanet (substellar, planetary-mass) candidates from astrometry,
  if stated; otherwise the number of DR4 orbital solutions whose catalogue or derived companion mass
  is below 13 MJup, counted by the rule in the release documentation. **RIGHT** if 3,300 <= N <= 11,700
  (forecast +- 2 sigma), **LOW** or **HIGH** otherwise. If neither count can be formed, L1 is
  **UNGRADABLE** and we say so.

Summary statistics: hit rate of the binary forecasts (RIGHT / graded), separately for "detectable" and
"not detectable" forecasts, and the two Brier scores. With 20 binary forecasts these are descriptive
only; no significance is claimed. Count, frozen: 20 binary forecasts (9 "detectable": P1-P7, P18,
B1; 11 "not detectable": P8-P17, B2), 2 probabilistic (W1, W2), 1 inclination test (B3) and 1 yield
forecast (L1): 24 in all.

## 5. What would and would not be claimed

- A forecast graded WRONG means the DR4 catalogue did not show what the paper forecast. It does not by
  itself mean the planet is not there: the DR4 pipeline may be more conservative than the forecast
  model, and companions with periods near or beyond the 66-month baseline are hard to solve. We will
  say this next to every WRONG.
- No forecast is re-interpreted after DR4. If a paper's wording allows two readings, the reading
  frozen in section 2 is used.
- If DR4's table structure makes a rule impossible to apply as written, the rule will be adapted in a
  dated amendment that changes only I/O (which table, which column), before any target is graded,
  and the change will be disclosed.

## 6. Timeline

Code and tests are committed before 2026-12-02. The job is filed as soon as the archive serves
gaiadr4 tables; grading follows the same day where compute allows.

## Credit

Forecasts: A. L. Wallace, R. A. Wittenmyer, J. Horner; L. Parc, J. P. Faria, F. Bouchy and the CORALIE
team; T. A. Baycroft, A. H. M. J. Triaud, J. Sahlmann; C. Lammers, J. N. Winn. Benedict & Harrison
(2017) for the HD 202206 inclinations. Data: ESA Gaia (DPAC), SIMBAD (CDS).

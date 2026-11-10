# Test fixtures and how they are regenerated

The files in this directory are bundled with the package and uploaded with every `grok publish`.
At runtime they are read from `System:AppData/NCA/<name>`. The PKNCA reference suite
(`src/tests/reference-suite-tests.ts`) and the other reference tests load them and check that the
DataFrame pipeline matches PKNCA within the package tolerances: 0.1 % on AUC and Cmax, 0.5 % on λz,
1 % on derived parameters. Reference values were produced by PKNCA 0.12.1.

Paths below are relative to the root of the `public` repository.

## Copies of the sci-comp NCA fixtures

These files are **byte-identical** copies of the sci-comp NCA reference suite:

```
libraries/sci-comp/src/nca/__tests__/datasets/01_theoph.csv      -> packages/NCA/files/theoph.csv
libraries/sci-comp/src/nca/__tests__/datasets/02_indometh.csv    -> packages/NCA/files/indometh.csv
libraries/sci-comp/src/nca/__tests__/datasets/04_iv_infusion.csv -> packages/NCA/files/04_iv_infusion.csv
libraries/sci-comp/src/nca/__tests__/fixtures/01_theoph.json     -> packages/NCA/files/theoph-reference.json
libraries/sci-comp/src/nca/__tests__/fixtures/02_indometh.json   -> packages/NCA/files/indometh-reference.json
```

Re-copy them after any change to the sci-comp suite:

```bash
SC=libraries/sci-comp/src/nca/__tests__
DST=packages/NCA/files
cp "$SC/datasets/01_theoph.csv"      "$DST/theoph.csv"
cp "$SC/datasets/02_indometh.csv"    "$DST/indometh.csv"
cp "$SC/datasets/04_iv_infusion.csv" "$DST/04_iv_infusion.csv"
cp "$SC/fixtures/01_theoph.json"     "$DST/theoph-reference.json"
cp "$SC/fixtures/02_indometh.json"   "$DST/indometh-reference.json"
```

When the reference JSON gained fields:

- **sci-comp 0.7.0:** the moment and lag parameters `aumclast`, `aumcinf_obs`, `mrt`, `vss`, `tlag`,
  `pct_aumcextrap`. The PKNCA call mapping is in sci-comp's own `REGEN.md`.
- **sci-comp 0.10.0:** PKNCA's `span.ratio`, merged into each profile's `provenance` as `span_ratio`.
  These JSON files are the only oracle for span ratio below 2 (theoph 1, 9, 10 and indometh 1); the
  rat CSVs have no such rows.

Do not mistake `config.min_span_ratio: 2` for an output. It is a PKNCA input option; the value per
profile is `provenance.span_ratio`.

`04_iv_infusion-reference.json` is **derived**, not copied. sci-comp's `04_iv_infusion.json` has no
per-profile `provenance` block, which the package's row assertion needs. The derived file has the
same `parameters` values plus a `provenance` stub in which `lambda_z_n_points` and `c0_extrapolated`
are null (an IV infusion has no back-extrapolated C0, and the λz point-count check is skipped). It is
the only fixture that exercises the IV-infusion −T_inf/2 correction of MRT and Vss.

## Rat and batch-screen datasets (PKNCA validation corpus)

Inputs and PKNCA references from the package's PKNCA validation corpus. The corpus and its R
generator scripts are maintained outside this repository.

| File | Contents |
|---|---|
| `03_rat_iv.csv`, `03_rat_po.csv`, `03_rat_sc.csv` | dataset 03, one route per file |
| `rat_iv_nca_datagrok_match.csv`, `rat_po_…`, `rat_sc_…` | PKNCA references for dataset 03 |
| `03_rat_multiroute.csv` (+ `_nonmem`) | dataset 03, all routes in one table (per-row Route; NONMEM event format) |
| `rat_multiroute_nca_datagrok_match.csv`, `rat_multiroute_bioavailability.csv` | references for the multi-route table and F% |
| `05_multicompound_batch_screen.csv` | dataset 05: 20 compounds × 3 rats, IV bolus, no body-weight column (per-kg CL / Vz path) |
| `multicompound_nca_per_subject.csv` | PKNCA reference for dataset 05, keyed on `SubjectID` (unique across compounds) |
| `06_rat_tk28d_nonmem.csv` | dataset 06: 28-day TK in NONMEM event format |
| `rat_tk28d_nca_per_cell.csv`, `rat_tk28d_accumulation_ratios.csv`, `rat_tk28d_exposure_margin.csv` | references for dataset 06 (per composite cell, Rac, exposure margin) |

The rat and batch references are CSV with PKNCA dot-named columns (`auclast`, `cmax`, `lambda.z`,
`aucinf.obs`, `cl.obs`, `vz.obs`, …) and are loaded with `loadReferenceCsv`, not by the JSON path.
Dataset 06 is asserted per cell in `reference-nonmem-tests.ts`, through the production
composite-pool path.

**BLQ policy.** The generators for datasets 05 and 06 pass
`conc.blq = list(first = "keep", middle = "keep", last = "keep")` explicitly, matching the package
default (BLQ set to zero before the first, between and after the last measurable sample). PKNCA
0.12.1's own default is `middle = "drop"`, so the option is pinned rather than left to PKNCA. On
these two datasets it changes nothing (measured: 0 differing cells), because every BLQ row is
trailing. The embedded-BLQ case is covered by a synthetic test in `reference-suite-tests.ts` and by
`06_blq_rules` in sci-comp.

## Sex-stratifier fixtures

Each is derived from a bundled source by rewriting **only the Sex column**; every other cell is
byte-identical (CRLF kept). An explicit value → M/F map collapses each variant back to the source's
strata, so the per-sex F% and TK reference values are reproduced exactly.

| Fixture | Source | Sex-column change |
|---|---|---|
| `03_rat_multiroute_gender.csv` | `03_rat_multiroute.csv` | header `Sex` → `Gender`, values unchanged |
| `03_rat_multiroute_sexcollision.csv` | `03_rat_multiroute.csv` | per subject: males in {`M`, `Male`, `MALE`}, females in {`F`, `Female`} |
| `06_rat_tk28d_nonmem_sexcoded.csv` | `06_rat_tk28d_nonmem.csv` | `M` → `1`, `F` → `2` |

The `1 = M, 2 = F` coding is ours, so the tests map it by construction; there is no universal
numeric sex code.

## Other references

- `theoph-exclusion-reference.json`: PKNCA 0.12.1 on R's built-in `Theoph` with the λz-dependent
  parameters of subjects 1 (span ratio < 2) and 5 (analyst reason) removed with PKNCA `exclude()`.
  MEAN / SD / N are over the rows `exclude()` leaves, as PKNCA's `summary()` aggregates them.
  Asserted by `src/tests/workflow-stages-exclusion-tests.ts`, which reads this file.
- `mouse_sparse-reference.json`, `*_nonmem-reference.json`, `04b_…-reference.json`: references for
  the sparse and NONMEM event-format tests.
- `07_iv_bolus_sparse.csv` + `07_iv_bolus_sparse-reference.json`: SIMULATED destructive IV-bolus
  composites (four dose groups) and their PKNCA 0.12.1 oracle (`pk.calc.c0`, `pk.calc.auciv`). Written,
  not copied, by sci-comp `src/nca/__tests__/regen-sparse-fixture.R` block 4 (seeded generator; run it
  from `libraries/sci-comp`); identical to sci-comp's `datasets/07_iv_bolus_sparse.csv` and
  `fixtures/07_iv_bolus_sparse.json`. Asserted by `src/tests/sparse-iv-anchor-tests.ts`.

## Demo tables (`demo/`)

Opened from the Hub's demo link. Each is a bundled source with **one appended `Smiles` column**; every
other cell, the line endings (CRLF) and the quoting are byte-identical to the source. They are not
validation fixtures.

| File | Source | Numbers |
|---|---|---|
| `demo/indomethacin_iv.csv` | `indometh.csv` + indomethacin `COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1` | **Real**: Kwan et al. 1976, R `datasets::Indometh` |
| `demo/propranolol_multiroute.csv` | `03_rat_multiroute.csv` + propranolol `CC(C)NCC(O)COc1cccc2ccccc12` | **Simulated**: one-compartment, seed 20260526, CL 0.5 L/h/kg, Vd 2.5 L/kg, F_PO 0.35, F_SC 0.85 |

The propranolol structure is **illustrative**. The simulated data are dose-linear and low-clearance;
real rat propranolol has saturable first pass and clearance of 43–94 mL/min/kg (Iwamoto et al. 1985,
DOI 10.1111/j.2042-7158.1985.tb03041.x). The demo discloses this in the demo balloon ("Simulated
data, illustrative structure."); the text and its four-word limit are pinned by a test. The frame tag
`nca.demo.provenance` marks both tables but is not shown anywhere in the platform UI.

Display names and the balloon text live in `DEMO_TABLE_FILES` (`src/views/nca-launch.ts`).

To check a regenerated demo table, strip the last field and diff against the source; the result must
be empty:

```bash
diff <(sed 's/,[^,]*$//' OLD.csv) <(sed 's/,[^,]*$//' packages/NCA/files/demo/NEW.csv)
```

## Shape notes

- **Theoph**: `Subject, Wt, Dose, Time, conc`. Dose is mg/kg (dose in mg = `Dose * Wt`); oral.
- **Indomethacin**: `Subject, time, conc`. Fixed 25 mg dose, given in the test metadata; IV bolus.
- **IV infusion** (`04_iv_infusion.csv`): `Subject, time, conc`. 100 mg over 1 h (scalars).
  Synthetic one-compartment model (CL 2, V 10): analytic AUCinf 50, MRT 5, Vss 10.
- **Reference JSON** keys are snake_case: `cmax, tmax, auclast, aucinf, pct_aucextrap, lambda_z,
  half_life, cl, vz`, plus `aumclast, aumcinf_obs, mrt, vss, tlag, pct_aumcextrap` (`vss` null for
  extravascular, `tlag` null for IV). Provenance scalars: `lambda_z_n_points, lambda_z_r_squared,
  lambda_z_adj_r_squared, lambda_z_time_first, lambda_z_time_last, clast_obs, c0_extrapolated`. There
  is only the point count, no list of points used.

## Demo studies (`demo-studies/`)

**All five tables are simulated.** The SMILES are real molecules with fictitious PK. NOAEL
(30 mg/kg/day), the projected human dose (150 mg) and body weight (70 kg), and so the exposure margin,
are design constants. `ncaSeedDemoStudies()` installs them as Hub studies (see `README.md`,
*Demo data → Demo studies*).

The tables come from a deterministic generator (Python standard library only) maintained outside this
repository. Accepted seeds: Study 1 `20261108`, Study 2 `20261104`, Study 3 `20261106`,
Study 4 `20261108`, Study 5 `20261104`.

| File | Study |
|---|---|
| `demo_study1_series_x_rat_iv_po.csv` | Series X rat IV/PO screen: 8 compounds × (IV 1 mg/kg + PO 5 mg/kg), n = 3, male SD |
| `demo_study2_dgk0007_mouse_sparse_po.csv` | DGK-0007 mouse PO 10 mg/kg, destructive, 3 / time / sex, actual times ±8 % |
| `demo_study3_dgk0007_rat_dose_escalation.csv` | DGK-0007 rat IV 0.3 / 1 / 3, PO 3 / 10 / 30 mg/kg, n = 3, male |
| `demo_study4_dgk0007_rat_tk28d.csv` | DGK-0007 28-day TK: PO 0 / 10 / 30 / 100 mg/kg/day, 3 / sex / dose, Day 1 + Day 28, τ = 24 h; subject id `<animal>_D<day>` |
| `demo_study5_dgk0007_rat_iv_infusion.csv` | DGK-0007 0.5 h IV infusion, 1 and 3 mg/kg, 3 / sex / dose |

**Model.** Two-compartment disposition (mouse: one-compartment) with first-order absorption, RK4.
Hepatic blood flow: rat 55, mouse 90 mL/min/kg; blood:plasma = 1. F = Fa·Fg·(1 − CL/Q_H) with
CL ≤ 0.73·Q_H (a between-animal draw is capped at 0.85·Q_H).

- **Lead DGK-0007:** CL 14 (male) / 9 (female) mL/min/kg, Vss 2.4 and Vz 3.4 L/kg, t½ 2.8 / 4.4 h.
- **Dose–exposure (Studies 3–5):** linear CL; solubility-limited absorption above 10 mg/kg
  (Fa·Fg 0.60 / 0.42 / 0.25 and ka 1.5 / 0.8 / 0.6 /h at ≤ 10 / 30 / 100 mg/kg).
- **TK nonlinearity:** time-dependent CYP inhibition at 100 mg/kg only (CL −35 % by Day 28), so
  Rac(AUCτ) ≈ 1.7 at 100 mg/kg and ≈ 1.0 at 10 / 30 mg/kg.
- **Distribution:** α = 2.5 /h for every compound.
- **Variability:** CL CV 20 %; V = f_CL^0.7 · e^η (η CV 10 %); inter-occasion CL CV 5 %. Day-28
  profiles carry the four preceding daily doses.
- **Residual error:** 10 % proportional + 0.2 ng/mL additive.
- **BLQ:** LLOQ 1 ng/mL; BLQ rows carry `BLQ = True` and `0.00`.
- **Vehicle controls (Study 4):** sampled at pre-dose, 2 h and 24 h on both days; all BLQ.

**Schedules.**

- IV: 0, 0.083, 0.25, 0.5, 1, 2, 4, 6, 8, 12, 24 h.
- PO: 0, 0.25, 0.5, 1, 2, 4, 6, 8, 12, 24 h.
- Infusion: 0, 0.25, 0.5 (end of infusion), 0.75, 1, 1.5, 2, 4, 6, 8, 12, 24 h.
- TK: 0, 0.5, 1, 2, 4, 6, 8, 12, 24 h.

The 12 h sample keeps a profile that is BLQ at 24 h from fitting λz over a span of less than two
half-lives.

**Planted review cases.** A planted animal has no between-animal variability and a 3 % residual, so
the planted feature decides its verdict.

| Ref | Profile(s) | Planted | Verdict | Pre-recorded |
|---|---|---|---|---|
| C1 | S1 DGK-0008, all 6 | t½ ≈ 15 h, last sample 24 h | Flagged span + %AUCextrap (20.9–29.8 %) | bulk Accept as is |
| C2 | S1 DGK-0001 PO, 3 | BLQ from 6 h | Flagged N ≤ 3 | Accept as is |
| C3 | S1 DGK-0006 PO R604 | second absorption input at 6 h | Not estimated (engine `partial`) | Acknowledge |
| C4′ | S4 vehicle controls, 12 profiles | 0 mg/kg, all BLQ | No exposure | Acknowledge (Expected (control)) |
| C5 | S3 PO 30 mg/kg R316 | alternating ×1.77 / ×0.23 on 6–24 h | Flagged adj. R² 0.886 | Exclude λz |
| C6 | S4 F 10 mg/kg Day 28 T010_D28 | 8 h result ×3.75 (lab-flagged, `Comment`) | Flagged adj. R² 0.890 | ✎ force-out 8 h → Accepted |
| C7 | S4 100 mg/kg Day 1 T019_D1, T022_D1 | absorption lag 3.5 h: Tmax 4 h (group 0.5–2 h) | Flagged N ≤ 3 (λz 8–24 h) | To review |
| C8 | S4 F 30 mg/kg Day 28 T016_D28 | alternating ×1.575 / ×0.425 on 6–24 h | Flagged adj. R² 0.893 | To review |
| C9 | S5 F 3 mg/kg R510 | 4 h ×1.04, 6 h ×0.97 | Flagged N ≤ 3 (auto-fit 8–24 h) | ✎ force-in 4 h, 6 h → 5 points, Accepted |
| C10 | S5 F 1 mg/kg R505 | end-of-infusion sample at 0.53 h | Accepted | — |
| C11 | S2 both composites | destructive design, 24 h all BLQ | Flagged N ≤ 3 | Accept as is |
| — | every other profile Flagged only by N ≤ 3 | — | Flagged N ≤ 3 | bulk Accept as is per analysis |

Notes:

- **N ≤ 3 is common on ordinary data.** The best fit keeps a 3-point window on 30–41 % of normal
  profiles at 5–10 % residual error (44 of 108 on this directory's rat and batch fixtures). Those
  profiles get one bulk Accept as is per analysis, and the seed asserts adj. R² ≥ 0.98 on each.
- **Vehicle controls in the TK analysis.** They are No exposure, with CL / Vz / `DOSE_NORM` empty,
  and no NaN or ±Infinity appears. Two known gaps: the Accumulation and Sex-ratio grids show an
  all-blank 0 mg/kg row, and the by-group table prints SD = 0.000 at N = 0.
- **Study 5 infusion** is stored in the recipe as a uniform 0.5 h value, because the analysis view
  ignores a per-row infusion-duration column. The `InfusionDuration_hr` column stays in the table,
  and the seed asserts that every row is 0.5 h.

**Screening.** The authoritative screen is the package test `Demo studies: planted verdicts`, run
on the shipped engine. It checks:

- C1–C11 as above;
- no verdict other than Accepted or N ≤ 3 only;
- every TK profile measurable at 24 h;
- span ≥ 2 on every dosed TK profile;
- Rac(AUCτ) in [1.4, 2.0] at 100 mg/kg and in [0.85, 1.15] at 10 / 30 mg/kg;
- Day-28 male t½ in [3.5, 5.5] h at 100 mg/kg;
- F ≤ 100 %;
- the Series X F order where adjacent compounds differ by more than 2 × SE.

**Series X structures** (IDs fictional, molecules real):

| ID | SMILES | CL (mL/min/kg) | Fa·Fg | F |
|---|---|---|---|---|
| DGK-0001 | `COc1ccc(cc1)C(=O)Nc1ccc(cc1)N1CCN(C)CC1` | 40 | 0.50 | 14 % |
| DGK-0002 | `O=C(Nc1ccc(cc1)N1CCN(C)CC1)c1ccccc1` | 30 | 0.60 | 27 % |
| DGK-0003 | `Clc1ccc(cc1)C(=O)Nc1ccc(cc1)N1CCN(C)CC1` | 22 | 0.60 | 36 % |
| DGK-0004 | `FC(F)(F)c1ccc(cc1)C(=O)Nc1ccc(cc1)N1CCN(C)CC1` | 16 | 0.60 | 43 % |
| DGK-0005 | `FC(F)(F)Oc1ccc(cc1)C(=O)Nc1ccc(cc1)N1CCN(C)CC1` | 12 | 0.30 | 23 % |
| DGK-0006 | `OC(=O)c1ccc(cc1)C(=O)Nc1ccc(cc1)N1CCN(C)CC1` | 6 | 0.06 | 5 % |
| DGK-0007 | `Fc1cc(NC(=O)c2ccc(cc2)C(F)(F)F)ccc1N1CCN(C)CC1` | 14 | 0.60 | 45 % |
| DGK-0008 | `CS(=O)(=O)c1ccc(cc1)C(=O)Nc1ccc(cc1)N1CCN(C)CC1` | 2.5 | 0.35 | 33 % |

The seed types the `SMILES` column as a molecule itself, so the demo installs without the Chem
package. Chem's detector types a SMILES column asynchronously, so a detect-then-check guard would
refuse every study.

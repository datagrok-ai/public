# Case 2 — Covariate Adjustment of Organ Weight by Body Weight (TOX-102, 3-month mouse)

> Demo dataset for the box-plot group-comparison tool: liver weight from NTP TOX-102,
> illustrating that organ weight cannot be compared as a raw mean — the conclusion depends
> on the normalization (absolute / relative / ANCOVA).

## 1. Source file (NTP CEBS)

Study: **TOX-102**, 3-month toxicity study of trans-resveratrol (CASRN 501-36-0), B6C3F1/N mice, gavage. Accession `002-02772-0009-0000-9`.

File used (individual-animal, true Excel `.xlsx`, sheet `Data`):

- [Individual Animal Organ Weight Data](https://cebs.niehs.nih.gov/cebs/get_file/accno/002-02772-0009-0000-9/file/2032304_Individual_Animal_Organ_Weight_Data.xlsx)

Publication root: [cebs.niehs.nih.gov/cebs/publication/TOX-102](https://cebs.niehs.nih.gov/cebs/publication/TOX-102) · Data DOI: [10.22427/NTP-DATA-TOX-102](https://doi.org/10.22427/NTP-DATA-TOX-102)

> **Why mouse, not rat.** TOX-102 has four toxicity studies; individual organ weights exist for the two 3-month arms (perinatal Wistar Han rat and B6C3F1/N mouse). The clean "covariate + adaptive-vs-adverse" signal is in the mouse: relative liver weight rises significantly in females from ≥625 mg/kg, and both absolute and relative liver weight rise in 2,500 mg/kg males, **with no microscopic lesions** (adaptive hypertrophy). Terminal body weight is already included in the organ-weight file, so the separate body-weight files are not needed.

## 2. Final CSV and how it was built

Final file: **[`tox102_3mo_mouse_organ_weights.csv`](tox102_3mo_mouse_organ_weights.csv)** (120 rows).

Steps (reproduced by [`build_tox102_case2.py`](build_tox102_case2.py)):

1. Download the organ-weight `.xlsx` (link above).
2. Read sheet `Data` with `pandas.read_excel(..., sheet_name='Data')`.
3. Rename columns, coerce weights and dose to numeric, recode `Male/Female → M/F`.
4. Compute relative liver weight `RelLiver_pct = Liver_g / TBW_g × 100`.
5. Sort, write CSV. No missing liver or body-weight values (N=120).

Run: `python3 build_tox102_case2.py` (deps: pandas, openpyxl, numpy, scipy, statsmodels, matplotlib).

## 3. CSV columns

| Column | Type | Description |
|---|---|---|
| `Animal` | int | Animal number |
| `Sex` | str | `M` / `F` |
| `Dose_mgkg` | int | Dose, mg/kg/day: 0, 156, 312, 625, 1250, 2500 |
| `is_control` | bool | `True` for 0 mg/kg |
| `TBW_g` | float | Terminal body weight, g (**covariate**) |
| `Liver_g` | float | Liver weight, g (absolute) |
| `RelLiver_pct` | float | Relative liver weight, % of body weight |
| `Heart_g`, `Kidney_R_g`, `Lungs_g`, `Thymus_g`, `Testis_R_g` | float | Other organs (for extension) |
| `Removal Day` | str | Removal day (all `SD92`) |

Structure: 6 doses × 2 sexes × N=10 = 120 animals.

## 4. Downstream processing & analysis pipeline

Demonstrates that the normalization choice changes the conclusion (knobs A→B→C→D + covariate):

1. **(A) Group** — box = `Dose_mgkg`.
2. **(B) Facet by sex**; **(C) matched control** within sex.
3. **Covariate — three displays of the same liver data:**
   - **Absolute** `Liver_g` (naive);
   - **Relative** `RelLiver_pct = Liver/BW` (industry adjustment);
   - **ANCOVA**: `Liver_g ~ Dose + TBW_g` — body weight as covariate (recommended method, STP/Sellers 2007).
4. **(D) Tests:** Dunnett (each dose vs control) for absolute and relative; for ANCOVA, an F-test of the dose effect after adjusting for body weight (comparison of models with and without dose).

## 5. Outputs

**Headline (females):** the same data yields **three different conclusions**.

| Method | Conclusion (females) | Significance |
|---|---|---|
| Absolute liver weight | "no effect" | all doses `ns` (max +10.5% at 625, p=0.17) |
| Relative liver weight | "clear dose effect" | significant from 625 (`*`), 1250 (`***`), 2500 (`***`) |
| **ANCOVA (body-weight adjusted)** | **effect is real, not a body-weight artifact** | dose effect p=0.00025; body-weight slope p≈3e-11 |

Why: in females, body weight **falls** with dose (29.1→26.5 g), so the absolute liver increase is masked while the relative weight inflates. ANCOVA adjudicates: the liver is genuinely heavier than body weight predicts (significant covariate slope). In males, body weight is stable, so all three methods agree (liver genuinely enlarged; 2,500 mg/kg significant).

**Biological meaning:** per the TOX-102 report, this liver-weight increase had **no microscopic correlate** → adaptive hepatocellular hypertrophy, not toxic injury. Demo takeaway: a "statistically significant organ-weight effect" ≠ "adverse," and the normalization must not be chosen silently — it changes the answer.

**Figures:**
- [`tox102_case2_three_normalizations.png`](tox102_case2_three_normalizations.png) — females, 3 panels: A absolute (`ns`), B relative (`***` from 625), C ANCOVA scatter (dose effect survives body-weight adjustment).
- [`tox102_case2_covariate_scatter.png`](tox102_case2_covariate_scatter.png) — liver vs body weight, both sexes, color = dose; high-dose points sit above the line (liver heavier than expected).

---

## Files in this folder
- [`tox102_3mo_mouse_organ_weights.csv`](tox102_3mo_mouse_organ_weights.csv) — cleaned dataset
- [`build_tox102_case2.py`](build_tox102_case2.py) — reproducible pipeline
- [`tox102_case2_three_normalizations.png`](tox102_case2_three_normalizations.png) — main teaching figure
- [`tox102_case2_covariate_scatter.png`](tox102_case2_covariate_scatter.png) — covariate scatter

## Data license
NTP/CEBS data are U.S. Government public-domain.

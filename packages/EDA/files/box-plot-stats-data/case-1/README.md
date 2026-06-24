# Case 1 — Sex-Stratified Baseline (TR-594, 3-month rat)

> Demo dataset for the box-plot group-comparison tool: terminal body weight from NTP TR-594,
> illustrating why the control must be matched per sex.

## 1. Source files (NTP CEBS)

Study: **TR-594**, 3-month inhalation study in Hsd:Sprague Dawley SD rats, p-chloro-α,α,α-trifluorotoluene (CASRN 98-56-6). Accession `002-03160-0004-0000-6`.

Files used (individual-animal, format: HTML table with `.xls` extension):

- [Male Body Weight](https://cebs.niehs.nih.gov/cebs/get_file/accno/002-03160-0004-0000-6/file/1047201_Male_Individual_Animal_Body_Weight_Data.xls)
- [Female Body Weight](https://cebs.niehs.nih.gov/cebs/get_file/accno/002-03160-0004-0000-6/file/1047201_Female_Individual_Animal_Body_Weight_Data.xls)

Publication root: [cebs.niehs.nih.gov/cebs/publication/TR-594](https://cebs.niehs.nih.gov/cebs/publication/TR-594) · Data DOI: [10.22427/NTP-DATA-TR-594](https://doi.org/10.22427/NTP-DATA-TR-594)

> ⚠️ **Note on organ-weight files.** The `..._Organ_Weight_Data.xls` files for this study **contain no organ weights** — only Animal Number, Dose Group, Control Code and Final Body Weight (a duplicate of terminal weight). So Case 1 uses **terminal body weight**, and the covariate Case 2 (organ ÷ body weight) must use a different dataset (resveratrol TOX-102).

## 2. Final CSV and how it was built

Final file: **[`tr594_3mo_rat_terminal_bw.csv`](tr594_3mo_rat_terminal_bw.csv)** (120 rows).

Steps (fully reproduced by [`build_tr594_case1.py`](build_tr594_case1.py)):

1. Download the two body-weight `.xls` files (links above).
2. Parse each with `pandas.read_html(path)[1]` — the data is the second table (first = metadata).
3. Map text groups to doses: `Vehicle Control→0, 125ppm→125, …, 2000ppm→2000`.
4. Take terminal weight from column `Week 14` (= necropsy-day body weight, ~day 93–94).
5. Concatenate males + females into long format, add `Sex`, `is_control`, sort.
6. Write CSV. One female (2000 ppm, animal 1109) was sacrificed moribund on day 45 → no terminal weight (`NaN`); the record is **kept** but excluded from terminal-endpoint analysis.

Run: `python3 build_tr594_case1.py` (deps: pandas, lxml, numpy, scipy, matplotlib).

## 3. CSV columns

| Column | Type | Description |
|---|---|---|
| `Animal` | int | Animal number (unique within sex) |
| `Sex` | str | `M` / `F` |
| `Dose_ppm` | int | Dose in ppm: 0, 125, 250, 500, 1000, 2000 |
| `Dose_label` | str | Original group label (`Vehicle Control`, `125ppm`, …) |
| `Days_on_study` | int | Days on study (terminal ≈ 93–94; early removals lower) |
| `Terminal_BW_g` | float | Terminal body weight, g (Week 14). `NaN` for early removals |
| `Removal_reason` | str | `Terminal Sacrifice` / `Moribund Sacrifice` |
| `is_control` | bool | `True` for the 0 ppm group |

Structure: 6 doses × 2 sexes × N=10 = 120 animals (terminal available for 119).

## 4. Downstream processing & analysis pipeline

Demonstrates naive vs structured comparison (knobs A→B→C→D):

1. **(A) Group** — box = `Dose_ppm`.
2. **(B) Facet by sex** — separate M and F panels; do **not** pool sexes.
3. **(C) Matched control** — each dose box compared to the control of its **own** sex.
4. **(D) Tests** (within each sex):
   - **Dunnett** — each dose vs its own control (multiplicity-adjusted);
   - **Jonckheere–Terpstra** — rank test for a monotonic dose trend.
5. **Counter-example (naive)** — same tests on the pooled M+F data with a single control, to show the failure.
6. **Display** — box + jittered points (N=10; points show the true spread).

## 5. Outputs

**Headline:** the dose effect is **sex-divergent**, and pooling sexes destroys it.

| | Control | Trend (Jonckheere) | Top dose vs control |
|---|---|---|---|
| ♂ Males | 393 g (SD 22) | **p = 0.029, ↓ decrease** | 2000 ppm: −17 g |
| ♀ Females | 240 g (SD 15) | **p < 0.0001, ↑ increase** | 2000 ppm: +54 g *** |
| Pooled M+F (naive) | 316 g (**SD 79**) | p = 0.33 (no trend) | all doses `ns` |

When pooled: (1) control SD inflates from ~15–22 to **79** g due to bimodality, (2) the opposite male/female trends cancel → naive conclusion "no effect." Sex stratification recovers both real effects. This is a direct illustration of the key point: what matters is not the test (knob D) but the comparison structure (knobs B + C).

**Figures:**
- [`tr594_case1_naive_vs_correct.png`](tr594_case1_naive_vs_correct.png) — 3 panels: A naive pool (bimodal control), B males matched, C females matched, with significance stars.
- [`tr594_case1_trends.png`](tr594_case1_trends.png) — dose response M vs F (± SE); dotted pooled line is nearly flat.

---

## Files in this folder
- [`tr594_3mo_rat_terminal_bw.csv`](tr594_3mo_rat_terminal_bw.csv) — cleaned dataset
- [`build_tr594_case1.py`](build_tr594_case1.py) — reproducible pipeline
- [`tr594_case1_naive_vs_correct.png`](tr594_case1_naive_vs_correct.png) — main teaching figure
- [`tr594_case1_trends.png`](tr594_case1_trends.png) — dose-response divergence

## Data license
NTP/CEBS data are U.S. Government public-domain.

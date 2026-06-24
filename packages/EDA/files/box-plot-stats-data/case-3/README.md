# Case 3 — Small N + Dose-Response Trend (CJ16050, SEND dataset, respiratory function in rats)

> Demo dataset for the box-plot group-comparison tool: respiratory endpoints from the PHUSE
> CJ16050 SEND example study. Illustrates two things: (1) at small N a strip plot is more honest
> than a box plot; (2) a trend test (Jonckheere) only sees monotonic effects — non-monotonic
> ones need pairwise Dunnett.

## 1. Source files (PHUSE, GitHub)

Study: **CJ16050** ("RE Function in Rats") — a CDISC SEND example dataset from the PHUSE phuse-scripts repository. Respiratory function, rats, Compound A, gavage.

SEND domains used (format `.xpt`, SAS Transport v5):

- `dm.xpt` (demographics: SEX, SETCD) — [raw.githubusercontent.com/…/dm.xpt](https://raw.githubusercontent.com/phuse-org/phuse-scripts/master/data/send/CJ16050/dm.xpt)
- `re.xpt` (respiratory endpoints: RESTRESN, RETESTCD, RETPT) — [raw.githubusercontent.com/…/re.xpt](https://raw.githubusercontent.com/phuse-org/phuse-scripts/master/data/send/CJ16050/re.xpt)

Full folder: [github.com/phuse-org/phuse-scripts/…/CJ16050](https://github.com/phuse-org/phuse-scripts/tree/master/data/send/CJ16050) (also contains ts/tx/ta/te/ex/se/cl/ds and define.xml — not needed here).

> **Study structure.** 3 dose groups (0 / 100 / 1000 mg/kg Compound A) × **6 males** = 18 rats. **Males only** — so sex faceting (knob B) is not applicable here; it is covered by Cases 1–2. This case's job is small N and dose trend. RE endpoints: Respiratory Rate (breaths/min), Tidal Volume (mL/breath), Minute Volume (mL/min) at Predose / 1 / 2 / 4 / 8 h postdose (Day 1).

## 2. Final CSV and how it was built

Final file: **[`cj16050_resp_function.csv`](cj16050_resp_function.csv)** (267 rows; 18 animals × 3 tests × 5 timepoints − missing).

Steps (reproduced by [`build_cj16050_case3.py`](build_cj16050_case3.py)):

1. Download `dm.xpt` and `re.xpt` (links above). GitHub is directly reachable via `curl`/`raw.githubusercontent.com`.
2. Read the domains with `pandas.read_sas(path, format='xport', encoding='latin1')` (XPORT is read natively, no extra library).
3. Join RE to DM on `USUBJID`; take dose from `SETCD` (`00→0, 01→100, 02→1000`; control is encoded as ≈0 in EX).
4. Coerce `RESTRESN` to numeric, reshape to long format (test × timepoint × animal).
5. Drop missing, write CSV.

Run: `python3 build_cj16050_case3.py` (deps: pandas, numpy, scipy, matplotlib).

## 3. CSV columns

| Column | Type | Description |
|---|---|---|
| `USUBJID` | str | Unique subject ID (SEND) |
| `Animal` | str | Short ID (last 3 chars of USUBJID) |
| `Sex` | str | `M` (all male) |
| `Dose_mgkg` | int | Dose, mg/kg: 0, 100, 1000 |
| `Dose_label` | str | Group label (`Control`, `Compound A 100 mg/kg`, …) |
| `TestCD` / `Test` | str | Test code/name: `RESPRATE` / `TIDALVOL` / `MV` |
| `Timepoint` / `TimepointNum` | str / float | Timepoint (`Predose`, `1-hour postdose`, …) and its number |
| `Value` | float | Numeric result (`RESTRESN`) |
| `Unit` | str | Unit (breaths/min, mL/breath, mL/min) |
| `is_control` | bool | `True` for 0 mg/kg |

## 4. Downstream processing & analysis pipeline

1. **(A) Group** — box = `Dose_mgkg` (3 ordered doses).
2. **(B) Facet by sex** — not applicable (males only).
3. **(C) Control** — the 0 mg/kg group (SETCD 00).
4. **Small N (N=6) — display choice:** with 6 points, box-plot quartiles are unreliable → show a **strip/scatter** (every animal) + mean/median.
5. **(D) Tests:** for one endpoint/timepoint (demo uses 4-hour postdose):
   - **Jonckheere–Terpstra** — rank test for a monotonic dose trend;
   - **Dunnett** — each dose vs control (pairwise).
   The right test depends on the **shape** of the response, which the strip plot makes visible.

## 5. Outputs

**Two lessons at 4-hour postdose:**

| Endpoint | Response shape | Jonckheere (trend) | Dunnett (pairwise) |
|---|---|---|---|
| Tidal Volume | monotonic decrease (1.63→1.55→1.00) | **p=0.003** ✓ detects | 100 `ns`, 1000 `***` |
| Respiratory Rate | biphasic (118→185→79) | **p=0.12** ✗ misses | 100 `***`, 1000 `*` |
| Minute Volume | biphasic (191→288→79) | p=0.075 ✗ weak | 100 `**`, 1000 `**` |

Demo takeaways:
1. **Small N:** a box plot on 6 points implies precision that is not there; the strip plot reveals the real spread — e.g., the 1000 mg/kg group has 1–2 "non-responders" (1.18–1.43 mL) against the rest at ~0.85, which the box hides.
2. **A trend test is not a silver bullet:** Jonckheere is sensitive only to a **monotonic** response (Tidal Volume) and misses the **biphasic** one (Respiratory Rate / Minute Volume), where Dunnett catches both the mid-dose increase and the high-dose decrease. The correct test choice is only visible after honestly displaying the response shape.

**Figures:**
- [`cj16050_case3_box_vs_strip.png`](cj16050_case3_box_vs_strip.png) — Tidal Volume @4h: box (false precision) vs strip (real spread), N=6.
- [`cj16050_case3_trend_vs_pairwise.png`](cj16050_case3_trend_vs_pairwise.png) — monotonic (Jonckheere detects) vs biphasic (Jonckheere misses, Dunnett catches).

---

## Files in this folder
- [`cj16050_resp_function.csv`](cj16050_resp_function.csv) — cleaned dataset
- [`build_cj16050_case3.py`](build_cj16050_case3.py) — reproducible pipeline
- [`cj16050_case3_box_vs_strip.png`](cj16050_case3_box_vs_strip.png) — small N: box vs strip
- [`cj16050_case3_trend_vs_pairwise.png`](cj16050_case3_trend_vs_pairwise.png) — monotonic vs biphasic trend

## Caveats and license
- **MIT License** (PHUSE phuse-scripts) — free for commercial/private use provided the license notice is retained. This is **synthetic/anonymized** example data, not a real regulatory submission.
- Doses derived from `SETCD` (control = 0; in `EX` control is encoded as ≈0, a near-zero artifact).
- N=6 per group; males only.

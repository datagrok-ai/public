# Reference fixture regeneration (PKNCA oracle)

The JSON files in `fixtures/` are the **binding reference oracle** for the NCA
reference suite (`core/__tests__/reference-suite.test.ts`). Every committed
parameter value is produced by **PKNCA 0.12.1** so the TypeScript core is
validated against an independent, peer-reviewed implementation (PRD NFR-04/05;
nca-studio CLAUDE.md rule 18). A parameter without a passing reference assertion
does not exist.

## Toolchain

- R ≥ 4.6.0
- `PKNCA` **0.12.1** (pinned — `regen-fixtures.R` asserts the version)
- `jsonlite`
- Node ≥ 18 (for the merge step)

## Procedure

```bash
# from libraries/sci-comp/
Rscript src/nca/__tests__/regen-fixtures.R     # 1. compute references
node    src/nca/__tests__/merge-fixtures.mjs    # 2. merge into 01/02/03
```

Step 1 (`regen-fixtures.R`):
- re-derives the **8 original parameters** (Cmax … Vz) and prints the max
  relative error vs the committed fixtures. This is a self-check: a non-zero
  error (> ~1e-9) means the PKNCA configuration here no longer reproduces the
  original run, so the **new** moment/lag values would not be trustworthy.
  Investigate before trusting the output. Current reproduction: **< 3e-11**.
- writes `fixtures/_new_params.json` — the 6 moment/lag fields **plus
  `span_ratio`** per subject (a transient intermediate; not itself an asset —
  do not commit it).
- writes `datasets/04_iv_infusion.csv` + `fixtures/04_iv_infusion.json` — the
  new IV-infusion fixture (see below).

Step 2 (`merge-fixtures.mjs`) injects `aumclast, aumcinf_obs, mrt, vss, tlag,
pct_aumcextrap` into each profile's `parameters`, and `span_ratio` into its
`provenance`, in `fixtures/0{1,2,3}.json` — **preserving every existing value
exactly** (it never recomputes the original 8).

### `span_ratio` (terminal-phase span, PKNCA `span.ratio`)

Lives in `provenance`, beside the other `lambda_z_*` fields, because it is a fit
diagnostic rather than a reported PK parameter. It is requested from PKNCA
directly (`span.ratio` is a real interval column) rather than derived from
`lambda_z_time_first/last` + `half_life`, so the fixture carries the number PKNCA
itself reports — the rule-18 oracle for `LambdaZResult.spanRatio`.

`reference-suite.test.ts` asserts every theoph and indometh subject against it,
and pins the corpus incidence the diagnostic exists for: **4 of 18 profiles sit
below the conventional 2**, worst `02_indometh` subject 1 at 0.685 with adjusted
R² 0.994.

## PKNCA configuration (matches each fixture's `config` block)

```r
PKNCA.options(
  auc.method = "lin up/log down",
  min.hl.points = 3,
  min.hl.r.squared = 0.85,
  allow.tmax.in.half.life = FALSE,   # exclude_cmax = TRUE
  min.span.ratio = 2
)
```

## Parameter mapping (PKNCA → core `ParameterValues`)

| Core field      | PKNCA param (EV)  | PKNCA param (IV)        |
|-----------------|-------------------|-------------------------|
| `aumcLast`      | `aumclast`        | `aumclast`              |
| `aumcInf`       | `aumcinf.obs`     | `aumcinf.obs`           |
| `mrt`           | `mrt.obs`         | `mrt.iv.obs` (−T_inf/2) |
| `vss`           | — (NaN, gated)    | `vss.iv.obs`            |
| `tlag`          | `tlag`            | — (null, gated)         |
| `pctExtrapAumc` | `(aumcinf.obs − aumclast)/aumcinf.obs·100` | same |

**Route handling per dataset** (PKNCA function per route, R1-F6):
- **Extravascular** (theoph, rat): `route = "extravascular"`. A predose
  `(t=0, conc=0)` row is inserted when absent — the same convention `computeNca`
  uses — so PKNCA integrates AUC/AUMC over `[0, Inf)`. `mrt.obs`,
  `duration.dose = 0`. `vss` is **not** requested (IV-only); `tlag` is taken
  from PKNCA.
- **IV bolus** (indometh): `route = "intravascular"`, `duration = 0`. The
  log-linear `c0` from the committed fixture provenance (`c0_extrapolated`,
  identical to the core's `insertC0`) is inserted at `t=0` so PKNCA integrates
  over the **same augmented profile** the core uses. `mrt.iv.obs`,
  `vss.iv.obs`. `tlag` is forced **null** — it is an absorption concept, N/A
  for IV (the core gates it to `NaN`; PKNCA's value on inserted-c0 data is a
  noise artifact).
- **IV infusion** (04): `route = "intravascular"`, `duration = T_inf = 1 h`.
  `mrt.iv.obs` carries the `−T_inf/2` correction; `vss.iv.obs = CL·MRT`.

## The IV-infusion fixture (04) — why synthetic

No existing fixture exercises the `−T_inf/2` correction, and validating an
infusion against an IV-bolus profile is impossible. `04_iv_infusion` is
simulated from a known **one-compartment** model so the truth is analytic:

```
CL = 2 L/h, V = 10 L  →  k = CL/V = 0.2 /h ;  dose = 100 mg, T_inf = 1 h
C(t≤T) = (R0/CL)(1 − e^{−kt}),  R0 = dose/T_inf
C(t>T) = C(T)·e^{−k(t−T)}
⇒ analytic  AUCinf = dose/CL = 50,  MRT_iv = 1/k = 5,  Vss = V = 10
```

PKNCA recovers AUCinf ≈ 49.99, MRT ≈ 5.003, Vss ≈ 10.008 — these PKNCA values
(not the analytic ones) are stored as the reference, so the suite tests
core-vs-PKNCA parity while the analytic values confirm the fixture itself is
sound.

## Validation output (measured, not assumed)

Max deviation of the TypeScript core (`computeNca`) vs the PKNCA references,
over **all four fixtures** (theoph, indometh, rat, infusion), in **both**
summation modes (naive and Neumaier-compensated):

| Parameter       | Max deviation        | Class metric        | Gate (AD-12) |
|-----------------|----------------------|---------------------|--------------|
| `aumcLast`      | 7.4e-12              | relative            | 1e-3         |
| `aumcInf`       | 1.9e-11              | relative            | 1e-3         |
| `mrt`           | 2.0e-11              | relative            | 1e-2         |
| `vss`           | 1.8e-11              | relative            | 1e-2         |
| `pctExtrapAumc` | 1.1e-9               | absolute (pp)       | 0.5 pp       |
| `tlag`          | 0 (exact)            | absolute            | exact        |

The core matches PKNCA to floating-point round-off — ~9 orders of magnitude
inside the dimensional-analogy gates — so the AD-12 tolerances hold with
enormous margin and were not widened. (Reproduce: temporarily compute the
per-parameter max deviation in `reference-suite.test.ts`; the gates encode the
class, the table above records the measurement.)

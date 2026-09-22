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
- writes `fixtures/06_blq_rules.json` — the per-rule BLQ oracle from
  `datasets/06_blq_rules.csv` (see "The BLQ-rules fixture"), and prints PKNCA's
  own `c0` for indometh plus the stock-PKNCA `auclast` contrast.

Step 2 (`merge-fixtures.mjs`) injects `aumclast, aumcinf_obs, mrt, vss, tlag,
pct_aumcextrap` into each profile's `parameters`, and `span_ratio` (all
datasets) plus `c0_pknca`, `pct_auc_back_extrap` (02 indometh only — a key the R
script did not produce for a dataset is left absent, never written as null)
into its `provenance`, in `fixtures/0{1,2,3}.json` — **preserving every
existing value exactly** (it never recomputes the original 8).

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

### `c0_pknca` and `pct_auc_back_extrap` (IV bolus only — 02 indometh)

`c0_pknca` is PKNCA's **own `c0` PPTESTCD** (`pk.calc.c0`, chain `c0 → logslope
→ c1 → cmin → set0`), requested on the **raw** indometh profiles (no inserted
row, `route = "intravascular"`, `duration = 0`). It is the independent oracle
for the core's back-extrapolation: the committed `c0_extrapolated` was the
core's own c0 fed back into PKNCA for the AUC run (see below), so until this
key existed nothing outside the core vouched for the number. Measured
2026-09-22: PKNCA's `c0` equals `c0_extrapolated` to ≥ 10 significant digits on
all 6 subjects, and it does NOT change when a `(0, 0)` pre-dose row is present
or when that row is numerically substituted (`conc.blq first = 0.01`) — PKNCA's
c0 ignores a substituted dose-time BLQ, which is the semantics `augmentProfile`
implements (a BLQ / non-positive dose-time row is *absent* for c0).

`pct_auc_back_extrap` is the back-extrapolated share of AUCinf — the Phoenix
WinNonlin `AUC_%Back_Ext_obs` analogue, which PKNCA does not report. It is a
stated formula on PKNCA inputs, not a PKNCA output: the dose-time → first-
observation segment `(c0 − C1)·t1 / ln(c0/C1)` (log-down, since `c0 > C1`;
linear `(c0 + C1)/2·t1` otherwise) with `c0 = c0_pknca`, over PKNCA's
`aucinf.obs`, × 100. Corpus: **every indometh subject back-extrapolates 16–28 %
of its AUCinf** (subject 1: 20.55 %). Diagnostic only; `computeNca` reports it as
`provenance.c0.pctAucBackExtrap`, asserted within `TOL.pctExtrap` (0.5 pp).

## The BLQ-rules fixture (06) — per-rule PKNCA oracle

`datasets/06_blq_rules.csv` (five hand-authored subjects, BLQ encoded as
`conc = 0` with `blq = 1`, single LLOQ 0.05) × `fixtures/06_blq_rules.json`
(FIVE rule blocks, each a real PKNCA run under the `conc.blq` option that means
the same thing as the sci-comp `BlqStrategy` it is paired with — the mapping is
in `regen-fixtures.R`):

| Block | sci-comp `BlqStrategy` | PKNCA `conc.blq` | Notes |
|---|---|---|---|
| R-A | `set-zero` × 4 | `keep` × 3 | the substituted 0 integrates; λz's own `conc ≤ 0` guard keeps it out of the fit |
| R-B | `exclude` × 4 | `drop` × 3 | the pre-GROK-20960 numbers — every rule used to integrate as `exclude` |
| R-C | `set-half-lloq` × 4 | `0.025` × 3 (numeric) | substitutes enter λz and move tlast |
| R-D | `missing` × 4 | BLQ rows physically removed (≡ `conc.na = "drop"`) | asserted identical to R-B on every parameter |
| R-E | nca-studio's shipped default `set-zero, set-zero, set-zero, exclude` | `keep` × 3 | PKNCA cannot express the afterLast split; the trailing run is trimmed either way, so R-E ≡ R-A — the GAP-W3 assertion |

Subjects: **P1** leading BLQ at 0 and 0.5 h, embedded BLQ at 6 h, positive tail;
**P2** leading BLQ, clean middle, two trailing BLQ; **P3** positive to 4 h then
three trailing BLQ with only two post-Cmax positives (λz NA under set-zero /
exclude / missing — the `status: 'partial'` oracle); **P4** clean log-linear
decay with one trailing BLQ that lands on the line under LLOQ/2; **I1** indometh
subject 1 + a FLAGGED `(0, 0, blq = 1)` pre-dose row (the dose-time gate under
every rule — identical numbers in all five blocks; c0 = PKNCA's own c0).

The regression pin: P1 `set-zero` AUClast 13.4729 vs `exclude` 15.8010 — the
"four rules bit-identical" defect cannot silently return.

### Measured divergences from PKNCA (asserted, not skipped)

FOUR cells in the 06 table are NOT parity, and `reference-suite.test.ts`
asserts what each side does from the fixture's own numbers:

1. **R-B AUC is NA in PKNCA.** `first = "drop"` removes the t=0 zero and PKNCA
   does not extrapolate to the interval start (the same NA as a raw IV-bolus
   profile), so `auclast` / `aucinf` are NA on every PO subject. sci-comp's
   `exclude` at t=0 drops the row and prepends `(0, 0)` by the extravascular
   convention — exactly the R-D construction, which is therefore the AUC oracle
   for `exclude` (`config.auc_oracle_block`). R-B pins Cmax / Tmax / λz / Tlag
   under PKNCA's own option.
2. **PKNCA reports λz fits below sci-comp's adj-R² floor.** `pk.calc.half.life`
   (0.12.1 source) selects by `lambda.z > 0` and the adj-R² tie-break only;
   `min.hl.r.squared` is not consulted there (it belongs to the post-hoc
   `exclude_half.life` helper). Under `set-half-lloq`, P2 and P3 get a PKNCA λz
   through the flat LLOQ/2 tail at r² 0.761 / 0.678 (adj 0.681 / 0.571);
   sci-comp's `minRSquared = 0.85` rejects every window → `'partial'`. A slope
   through substitutes is not a terminal phase; "not estimated" is the honest
   reading. The suite asserts the rejection AND that PKNCA's adj-R² is below
   the floor, and that these are the only two such cells.
3. **PKNCA's AUCinf under numeric substitution is two-profiled.** With trailing
   LLOQ/2 substitutes PKNCA integrates `auclast` through them (P4: to 24 h) but
   keeps `tlast` / `clast.obs` at the last above-LOQ observation (12 h, 0.32)
   and extrapolates `aucinf.obs = auclast + clast.obs / λz` from THAT. sci-comp
   honours the substitutes consistently (D1): `cLast = LLOQ/2` at the last
   sample, `AUCinf = AUClast + cLast/λz`. **Provenance of this choice: `house`,
   NOT `winnonlin`.** An earlier draft called it "the Phoenix behaviour"; peer
   review (2026-09-22) could not verify that against a primary Certara source,
   so the attribution is withdrawn rather than left asserted. The choice stands
   on its own merit — ONE terminal anchor serves both limbs (the `cLast`/`tLast`
   that ends AUClast starts the extrapolated tail), where PKNCA's treatment is
   internally mixed. Cmax, AUClast, AUMClast and the λz window still agree; the
   AUCinf family (AUCinf, CL, Vz, %extrap, AUMCinf, MRT) differs and is asserted
   on both sides.
4. **A substituted point can enter an ACCEPTED λz fit, and adjusted R² cannot
   say so.** P4 under R-C: the trailing LLOQ/2 substitute at t = 24 is positive,
   survives the trailing trim, joins the window and scores adj-R² 0.9995 — high
   *because* it sits near the line, not because it was measured. It also becomes
   the terminal anchor of the AUCinf tail. The fit statistics are structurally
   incapable of distinguishing this from a genuine observation, so the engine
   REPORTS it: `LAMBDAZ_SUBSTITUTED_BLQ` (severity `warning`) names the offending
   times whenever a point in the accepted window is BLQ-flagged but not dropped.
   Found by peer review, 2026-09-22; pinned by the P4 assertion in
   `reference-suite.test.ts` and four unit cases in `compute-nca.test.ts`.

Tlag: PKNCA computes `tlag` on the PRE-clean data, so it agrees with sci-comp's
mask-based definition under keep / drop / numeric (P1 = 0.5 h in R-A, R-B, R-C,
R-E). Only R-D (BLQ rows physically removed, predose 0 prepended) cannot see the
lag and reports 0; sci-comp stays mask-based (0.5 h) under every rule (AD-7).

### Exposure note: the route-aware gate widens how often `logslope` runs

Worth stating plainly, because it is a change in *population* rather than in
formula. Before GROK-20960 an IV-bolus profile carrying a non-positive pre-dose
row never reached the c0 chain at all — `hasT0` was true and the profile
integrated from the `(0, 0)`. Now such a profile is routed through
`c0 → logslope → c1 → cmin → set0`, and `logslope` extrapolates from the **first
two** post-dose measurable points only (`c0.ts`). For a genuinely biexponential
IV bolus whose first two samples straddle the distribution phase, a two-point
log-linear back-extrapolation can misestimate C0 — and with it CL, Vz and the
back-extrapolated share. The chain itself is unchanged (it is PKNCA's, ported),
but it is now *hit far more often*, and the committed fixtures exercise one
decay geometry (indometh subject 1), not a fast-distribution shape.

Two things bound the risk rather than remove it: `provenance.c0.pctAucBackExtrap`
reports how much of AUC₀–∞ actually rides on the extrapolation (16–28 % across
indometh), and `C0_FALLBACK` fires when the log-slope was not estimable at all.
Neither detects a *confidently wrong* two-point slope. Raised by peer review
2026-09-22; a diagnostic for widely-spaced early samples relative to the apparent
terminal half-life is filed as a follow-up, not shipped here.

### Provenance of the vendor attributions in this file (audited 2026-09-22)

Two peer-review rounds checked the "Phoenix / WinNonlin" claims here against
Certara's published documentation. The outcome was NOT uniform, so it is recorded
per claim rather than left to the reader to assume a project-wide audit:

| Claim | Status |
|---|---|
| `AUC_%Back_Ext` (`_obs`) is a real, named Phoenix parameter | **VERIFIED** (round 1, Certara Phoenix online help — NCA parameter formulas / discrete plasma parameters) |
| Phoenix integrates IV-bolus AUC FROM a back-extrapolated `C0`, with a 2-point log-linear method falling back to the first observed value | **VERIFIED** (round 1, same source) — this is the basis of the convention section below |
| A trailing LLOQ/2 SUBSTITUTE is used as `Clast`/`Tlast` for the AUCinf term | **WITHDRAWN** — could not be verified against a primary Certara source; re-tagged `house` and justified on self-consistency (divergence #3 above) |

The distinction is deliberate: a citation that was checked and held is kept, one
that could not be checked is withdrawn. Anything added here later gets the same
treatment — cite the section, or tag it `house`.

## IV-bolus AUC convention — sci-comp (WinNonlin) vs stock PKNCA

sci-comp integrates an IV-bolus profile **from the back-extrapolated c0**: the
augmented profile `(0, c0), (t1, C1), …` is what `computeNca` integrates and
fits λz on, so AUClast includes the dose-time → first-sample segment. That is
the Phoenix WinNonlin convention (`AUC` with `C0`, reported alongside
`AUC_%Back_Ext`). **Stock PKNCA does not do this.** Measured 2026-09-22
(PKNCA 0.12.1, `regen-fixtures.R` prints it on every run):

| indometh subject 1 | `auclast` |
|---|---|
| raw profile (no t=0 datum) | **NA** |
| with an observed `(0, 0)` row | **1.719365** (integrates from 0) |
| with `(0, 0)` + `conc.blq first = "drop"` | **NA** |
| sci-comp / committed fixture (from c0 = 2.3936) | **2.009898** |

This is why the 02 fixture feeds PKNCA the **augmented** profile (the core's c0
inserted at `t = 0`): PKNCA validates the *integration of that profile*, not the
*choice of convention*. The convention is a deliberate, documented design
decision (nca-studio `nca-calculation-standards.md` §7), not a fixture-validated
one, and must not be read as such. Consistency is the reason: the same subject
must not lose 14.5 % of its AUC because a pre-dose sample happened to be drawn —
before GROK-20960 a `(0, 0)` row was taken as a measured t=0 value and
integrated from 0 (`1.7194` above); `augmentProfile` now REPLACES a BLQ /
non-positive dose-time row by `(0, c0)`, so the with-row and no-row profiles
give identical parameters (asserted in `reference-suite.test.ts`, Fx-1).

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

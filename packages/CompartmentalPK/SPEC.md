# Application Specification: Multi-Compartment Pharmacokinetic Model with Brute-Force Parameter Optimization

<!-- TOC for partial reading -->
<!-- Section 1: General Architecture — lines 15–490 -->
<!-- Section 2: Main View — lines 491–500 -->
<!-- Section 3: Controls — lines 501–720 -->
<!-- Section 4: Display Elements — lines 721–840 -->
<!-- Section 5: Layout — lines 841–940 -->
<!-- Section 6: User Feedback — lines 941–1030 -->
<!-- Section 7: Validation — lines 1031–1110 -->
<!-- Section 8: Pipeline — lines 1111–1260 -->
<!-- Section 9: Reactivity — lines 1261–1330 -->
<!-- Sections 10–15: Data, Errors, Resources, Closure, UX, Testing — lines 1331–1620 -->

## 1. General Architecture

### 1.0. General Information

| Field | Value |
|---|---|
| Application name | Multi-Compartment PK Model |
| Package | CompartmentalPK |
| Entry function | `compartmentalPkApp()` |
| Brief description | Interactive two- and three-compartment pharmacokinetic model with IV bolus, IV infusion, and oral absorption input modes. Includes brute-force grid search parameter optimization against synthetic noisy observations, optimization landscape visualization, and animated compartment diagrams. |
| Main view type | `DG.TableView` |

### 1.1. Core

#### Task List

| Task ID | Name | Pipeline type | Trigger | Synchronicity | Execution environment | Parallelization |
|---|---|---|---|---|---|---|
| `task_primary` | PK ODE solution | Primary (reactive) | Any model/dosing/display input change | Sync | Main thread | No |
| `task_generate_data` | Generate synthetic observations | Secondary (on demand) | Button `btn_generate_data` | Sync | Main thread | No |
| `task_optimize` | Brute-force grid search | Secondary (on demand) | Button `btn_run_optimize` | Async | Web workers (parallel) | Yes — grid points distributed across worker pool of size `Math.max(1, navigator.hardwareConcurrency - 2)` |

#### Dependencies Between Tasks

```
task_primary        — independent; produces concentration-time solution
task_generate_data  — depends on task_primary (interpolates from primary solution at sampling time points, adds noise)
task_optimize       — depends on task_generate_data (needs observed data);
                      after completion, can feed best-fit params back to primary controls
```

---

#### Task: `task_primary`

**General characteristics:**

| Field | Value |
|---|---|
| Name | PK ODE solution |
| Description | Solves the two- or three-compartment PK ODE system with the selected administration route (IV bolus, IV infusion, or oral absorption), computing concentration-time profiles, compartment amounts, and derived PK metrics. Handles multiple dosing by segmented integration. |
| Dependency on other tasks | No |

**Computation Formulas and Model**

**Level 1 — required minimum:**

Variables:

| Variable | Meaning | Units | Domain |
|---|---|---|---|
| `A_c` | Drug amount in central compartment | mg | `A_c ≥ 0` |
| `A_p` | Drug amount in peripheral compartment | mg | `A_p ≥ 0` |
| `A_d` | Drug amount in deep compartment (3-comp only) | mg | `A_d ≥ 0` |
| `A_gut` | Drug amount in GI tract (oral only) | mg | `A_gut ≥ 0` |
| `CL` | Clearance | L/h | `CL > 0` |
| `V_c` | Central volume of distribution | L | `V_c > 0` |
| `V_p` | Peripheral volume of distribution | L | `V_p > 0` |
| `Q` | Intercompartmental clearance | L/h | `Q > 0` |
| `V_d` | Deep compartment volume (3-comp) | L | `V_d > 0` |
| `Q2` | Deep intercompartmental clearance (3-comp) | L/h | `Q2 > 0` |
| `k_a` | Absorption rate constant (oral) | h⁻¹ | `k_a > 0` |
| `F` | Bioavailability (oral) | dimensionless | `0 < F ≤ 1` |
| `Dose` | Administered dose | mg | `Dose > 0` |
| `tau` | Dosing interval | h | `tau > 0` |
| `n_doses` | Number of doses | integer | `n_doses ≥ 1` |
| `T_inf` | Infusion duration | h | `T_inf > 0` |
| `MEC` | Minimum effective concentration | mg/L | `MEC ≥ 0` |
| `MTC` | Minimum toxic concentration | mg/L | `MTC > MEC` |

Micro-constants derived from physiological parameters:

```
k₁₀ = CL / V_c                          (elimination rate)
k₁₂ = Q / V_c,    k₂₁ = Q / V_p        (central ↔ peripheral)
k₁₃ = Q₂ / V_c,   k₃₁ = Q₂ / V_d      (central ↔ deep, 3-comp only)
```

Two-compartment ODE system:

```
dA_c/dt = −k₁₀·A_c − k₁₂·A_c + k₂₁·A_p + R_in(t)
dA_p/dt =  k₁₂·A_c − k₂₁·A_p
```

Three-compartment ODE system:

```
dA_c/dt = −k₁₀·A_c − k₁₂·A_c − k₁₃·A_c + k₂₁·A_p + k₃₁·A_d + R_in(t)
dA_p/dt =  k₁₂·A_c − k₂₁·A_p
dA_d/dt =  k₁₃·A_c − k₃₁·A_d
```

Plasma concentration:

```
C(t) = A_c(t) / V_c
```

Input function R_in(t) by administration route:

| Route | R_in(t) | Initial conditions |
|---|---|---|
| IV Bolus | `R_in = 0` (dose applied as impulse: `A_c += Dose` at each dose time) | `A_c(0) = Dose`, `A_p(0) = 0`, `[A_d(0) = 0]` |
| IV Infusion | `R_in = Dose / T_inf` during infusion, `0` otherwise | `A_c(0) = 0`, `A_p(0) = 0`, `[A_d(0) = 0]` |
| Oral | `R_in = k_a · A_gut` with `dA_gut/dt = −k_a · A_gut` | `A_gut(0) = F · Dose`, `A_c(0) = 0`, `A_p(0) = 0`, `[A_d(0) = 0]` |

Multiple dosing (segment-by-segment integration):

```
For n_doses > 1:
  Integration is split at dose times: t_dose = [0, tau, 2·tau, ..., (n_doses−1)·tau]

  For each segment [t_dose_k, t_dose_{k+1}] (or [t_dose_last, T_end]):
    1. Apply dose event at t_dose_k:
       - IV Bolus: A_c += Dose
       - IV Infusion: set R_in = Dose/T_inf (active for T_inf hours)
       - Oral: A_gut += F · Dose
    2. Integrate ODE from t_dose_k to next event time
    3. Use end state as initial condition for next segment

  For IV Infusion, additional events at t_dose_k + T_inf (infusion end):
    - Set R_in = 0 for central compartment contribution from this dose
```

Derived PK metrics:

```
C_max  = max(C(t))
t_max  = t at which C(t) = C_max
AUC    = trapezoidal integration of C(t) over [0, T_end]

For multiple dosing (n_doses > 1):
  C_ss_max = max(C(t)) in the last dosing interval
  C_ss_min = min(C(t)) in the last dosing interval
  t_90_ss  = earliest t where C_max_interval / C_ss_max ≥ 0.9
             (time to reach ~90% of steady state)
```

Output properties (invariants):

| Property | Description |
|---|---|
| `A_c(t) ≥ 0` for all `t` | Drug amounts are non-negative |
| `A_p(t) ≥ 0` for all `t` | Peripheral amount is non-negative |
| `A_d(t) ≥ 0` for all `t` (3-comp) | Deep compartment amount is non-negative |
| `C(0) = Dose / V_c` for IV Bolus | Initial concentration from bolus |
| `C(0) = 0` for Infusion and Oral | No instantaneous concentration for non-bolus routes |
| Mass balance: total amount ≤ cumulative dose | Drug cannot be created |
| `C(t) → 0` as `t → ∞` | All drug is eventually eliminated |
| Multiple dosing: `C_ss_max > C_ss_min` | Steady-state trough is below peak |

Reference examples:

| # | Model | Route | Inputs | Expected output | Source |
|---|---|---|---|---|---|
| 1 | 2-comp | IV Bolus | CL=5, V_c=20, V_p=40, Q=10, Dose=100, t=0 | C(0)=5 mg/L, dA_c/dt=−75 mg/h, dA_p/dt=50 mg/h | Manual: k₁₀=0.25, k₁₂=0.5, k₂₁=0.25; dA_c/dt=−0.25·100−0.5·100+0.25·0=−75 |
| 2 | 3-comp | IV Bolus | CL=5, V_c=20, V_p=40, Q=10, V_d=60, Q₂=2, Dose=100, t=0 | C(0)=5 mg/L, dA_c/dt=−85 mg/h, dA_p/dt=50 mg/h, dA_d/dt=10 mg/h | Manual: k₁₃=0.1, k₃₁=0.033; dA_c/dt=−0.25·100−0.5·100−0.1·100=−85 |
| 3 | 2-comp | Oral | CL=5, V_c=20, V_p=40, Q=10, k_a=1.0, F=0.8, Dose=100, t=0 | A_gut(0)=80, dA_gut/dt=−80, dA_c/dt=80, dA_p/dt=0 | Manual: A_gut=0.8·100=80; R_in=1.0·80=80; dA_c/dt=−0.25·0−0.5·0+0.25·0+80=80 |
| 4 | 2-comp | IV Infusion | CL=5, V_c=20, V_p=40, Q=10, Dose=100, T_inf=1, t=0.5 | R_in=100 mg/h, A_c(0)=0, dA_c/dt|₀=100 mg/h | Manual: R_in=100/1=100; dA_c/dt=−0.25·0−0.5·0+0.25·0+100=100 |

**Level 2 — full formalization:**

- Complete mathematical formulation: multi-compartment linear mammillary PK model. The system is linear in state variables (given rate constants). For IV bolus and oral, dose events are handled as impulses (state variable increments) between integration segments. For infusion, R_in is a piecewise-constant forcing function.
- Analytical properties: for single-dose IV bolus, the two-compartment model has an analytical solution as a bi-exponential: `C(t) = A·exp(−α·t) + B·exp(−β·t)`, where α and β are eigenvalues of the rate constant matrix. This provides additional verification targets. The system is asymptotically stable (all eigenvalues negative).
- Numerical method: MRT (Modified Rosenbrock Triple) from `diff-grok` library. Suitable for stiff and non-stiff systems. Adaptive step size ensures accuracy at discontinuities (dose events) when integration is segmented appropriately.

Level 2 document location: in this specification.

**Task input parameters:**

| Parameter | Type | Units | Domain | Description |
|---|---|---|---|---|
| `modelType` | `'2-compartment' \| '3-compartment'` | — | — | Compartment model selection |
| `inputMode` | `'iv-bolus' \| 'iv-infusion' \| 'oral'` | — | — | Administration route |
| `cl` | `number` | L/h | `> 0` | Clearance |
| `vc` | `number` | L | `> 0` | Central volume |
| `vp` | `number` | L | `> 0` | Peripheral volume |
| `q` | `number` | L/h | `> 0` | Intercompartmental clearance |
| `vd` | `number` | L | `> 0` | Deep compartment volume (3-comp) |
| `q2` | `number` | L/h | `> 0` | Deep intercompartmental clearance (3-comp) |
| `ka` | `number` | h⁻¹ | `> 0` | Absorption rate (oral) |
| `f` | `number` | — | `(0, 1]` | Bioavailability (oral) |
| `dose` | `number` | mg | `> 0` | Dose amount |
| `repeatedDosing` | `boolean` | — | — | Whether to use multiple doses |
| `tau` | `number` | h | `> 0` | Dosing interval |
| `nDoses` | `number` | — | `≥ 1, integer` | Number of doses |
| `tInf` | `number` | h | `> 0` | Infusion duration |
| `tEnd` | `number` | h | `> 0` | Simulation end time |
| `mec` | `number` | mg/L | `≥ 0` | Minimum effective concentration |
| `mtc` | `number` | mg/L | `> 0` | Minimum toxic concentration |

**Task output data:**

| Parameter | Type | Description |
|---|---|---|
| `t` | `Float64Array` | Time values (h) |
| `ac` | `Float64Array` | Central compartment amount (mg) |
| `ap` | `Float64Array` | Peripheral compartment amount (mg) |
| `ad` | `Float64Array \| null` | Deep compartment amount (mg), null for 2-comp |
| `agut` | `Float64Array \| null` | Gut amount (mg), null for non-oral |
| `conc` | `Float64Array` | Plasma concentration C(t) = A_c(t)/V_c (mg/L) |
| `cMax` | `number` | Maximum concentration |
| `tMax` | `number` | Time of maximum concentration |
| `auc` | `number` | Area under concentration-time curve |
| `cssMax` | `number \| null` | Steady-state max concentration (null if single dose) |
| `cssMin` | `number \| null` | Steady-state min concentration (null if single dose) |
| `t90ss` | `number \| null` | Time to 90% steady state (null if single dose) |

**Computation implementation:**

| Step | Implementation method | Details | Documentation |
|---|---|---|---|
| 1 | Custom method | Compute micro-constants (k₁₀, k₁₂, k₂₁, k₁₃, k₃₁) from physiological parameters | Formulas in this spec |
| 2 | Custom method | Build segment schedule: list of integration intervals based on dose times and infusion events | Algorithm described in this spec |
| 3 | External library | `diff-grok` v1.2.0+, function `mrt(task: ODEs)`. Solve ODE for each segment. MRT is adaptive multi-rate, suitable for stiff/non-stiff. | [diff-grok](https://github.com/datagrok-ai/diff-grok) |
| 4 | Custom method | Concatenate segment solutions, compute derived metrics (C_max, t_max, AUC, steady-state metrics) | Formulas in this spec |

**Execution environment constraint:** main thread (synchronous, typically < 100 ms for standard parameter ranges).

---

#### Task: `task_generate_data`

**General characteristics:**

| Field | Value |
|---|---|
| Name | Generate synthetic observations |
| Description | Samples the primary ODE solution at predefined clinical time points and adds log-normal noise to simulate realistic PK observations. |
| Dependency on other tasks | Uses results of `task_primary` (concentration-time arrays for interpolation) |

**Computation Formulas and Model**

**Level 1 — required minimum:**

Variables:

| Variable | Meaning | Units | Domain |
|---|---|---|---|
| `C_true,i` | True concentration at sampling time `t_i` | mg/L | `C_true,i ≥ 0` |
| `C_obs,i` | Observed (noisy) concentration at `t_i` | mg/L | `C_obs,i > 0` |
| `CV` | Coefficient of variation for noise | % | `5 ≤ CV ≤ 40` |
| `σ` | Standard deviation of log-normal noise | dimensionless | `σ = sqrt(ln(1 + (CV/100)²))` |

Sampling schedules (time points in hours):

```
Standard PK: [0.5, 1, 2, 4, 8, 12, 24]
Rich:        [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 10, 12, 16, 20, 24]
Sparse:      [1, 4, 12, 24]
```

For multiple dosing: sampling times are extended to cover the full dosing schedule. Time points within each dosing interval are repeated relative to the last dose time.

Noise model (log-normal):

```
σ = sqrt(ln(1 + (CV/100)²))
C_obs,i = C_true,i · exp(N(0, σ²) − σ²/2)
```

Where `N(0, σ²)` is a normal random variable with mean 0 and variance σ². The `−σ²/2` correction ensures `E[C_obs] = C_true` (unbiased on the original scale).

Output properties:

| Property | Description |
|---|---|
| `C_obs,i > 0` | Observed concentrations are always positive (log-normal guarantee) |
| `E[C_obs,i] = C_true,i` | Noise is unbiased on original scale |

Reference examples:

| # | Inputs | Expected output | Source |
|---|---|---|---|
| 1 | CV=0%, C_true=5.0 | C_obs=5.0 (no noise) | By definition |
| 2 | CV=15%, many samples | mean(C_obs/C_true) ≈ 1.0 | Statistical property of log-normal |

**Task input parameters:**

| Parameter | Type | Units | Domain | Description |
|---|---|---|---|---|
| `primarySolution` | `PrimarySolution` | — | — | Solution from task_primary (t, conc arrays) |
| `cv` | `number` | % | `[5, 40]` | Coefficient of variation |
| `schedule` | `'standard' \| 'rich' \| 'sparse'` | — | — | Sampling schedule |
| `tEnd` | `number` | h | `> 0` | Simulation end time (for clipping sample times) |

**Task output data:**

| Parameter | Type | Description |
|---|---|---|
| `tObs` | `number[]` | Sampling time points |
| `cObs` | `number[]` | Observed (noisy) concentrations |
| `cTrue` | `number[]` | True concentrations at sampling times |

**Computation implementation:**

| Step | Implementation method | Details | Documentation |
|---|---|---|---|
| 1 | Custom method | Select sampling time points based on schedule; clip to T_end | — |
| 2 | Custom method | Interpolate primary solution at sampling times (linear interpolation between nearest ODE output points) | — |
| 3 | Custom method | Add log-normal noise using `Math.random()` with Box-Muller transform | — |

**Execution environment constraint:** main thread (synchronous, < 10 ms).

---

#### Task: `task_optimize`

**General characteristics:**

| Field | Value |
|---|---|
| Name | Brute-force grid search optimization |
| Description | Performs exhaustive grid search over 1–3 selected PK parameters. For each grid point, solves the full PK ODE and computes the WSSR against observed data. Returns the best-fit parameters, goodness-of-fit metrics, and the full WSSR landscape for visualization. |
| Dependency on other tasks | Uses observed data from `task_generate_data`; after completion, optionally feeds best-fit parameters back to primary controls |

**Computation Formulas and Model**

**Level 1 — required minimum:**

Variables:

| Variable | Meaning | Units | Domain |
|---|---|---|---|
| `C_obs,i` | Observed concentration at time `t_i` | mg/L | `> 0` |
| `C_pred,i` | Predicted concentration at time `t_i` | mg/L | `≥ 0` |
| `w_i` | Weight for observation `i` | varies | `w_i > 0` |
| `WSSR` | Weighted sum of squared residuals | (mg/L)² | `WSSR ≥ 0` |
| `R²` | Coefficient of determination | dimensionless | `R² ≤ 1` |
| `AIC` | Akaike Information Criterion | dimensionless | — |

Objective function:

```
WSSR = Σᵢ wᵢ · (C_obs,i − C_pred,i)²
```

Weight options:

```
Uniform:    wᵢ = 1
1/C_obs²:   wᵢ = 1 / C_obs,i²
1/C_pred²:  wᵢ = 1 / C_pred,i²
```

Goodness-of-fit metrics:

```
SS_tot = Σᵢ wᵢ · (C_obs,i − mean(C_obs))²
R² = 1 − WSSR / SS_tot

n = number of observations
p = number of optimized parameters
AIC = n · ln(WSSR / n) + 2·p
```

Grid search:

```
For each selected parameter p_k (k = 1..K, K ∈ {1,2,3}):
  grid_k = linspace(p_k_min, p_k_max, grid_resolution)

Total grid points = grid_resolution^K
For each combination (g_1, g_2, ..., g_K):
  Set selected parameters to grid values, keep others fixed
  Solve ODE → compute C_pred at observation times
  Compute WSSR

Best-fit = combination with minimum WSSR
```

Output properties:

| Property | Description |
|---|---|
| `WSSR_best ≤ WSSR` for all grid points | Best-fit minimizes objective |
| `R² ≤ 1` | By definition |
| Grid is exhaustive | Every combination is evaluated |

Reference examples:

| # | Inputs | Expected output | Source |
|---|---|---|---|
| 1 | Observed data generated with known params, 0% noise, optimize same params | Best-fit recovers original params exactly | By construction |
| 2 | 1 parameter, 20 grid points | WSSR landscape has single minimum | For identifiable parameter |

**Task input parameters:**

| Parameter | Type | Units | Domain | Description |
|---|---|---|---|---|
| `observedData` | `ObservedData` | — | — | `{tObs, cObs}` from task_generate_data |
| `selectedParams` | `SelectedParam[]` | — | 1–3 items | Parameters to optimize with min/max ranges |
| `fixedParams` | `PkParams` | — | — | All PK parameters not being optimized |
| `gridResolution` | `number` | — | `[5, 50]` | Grid points per axis |
| `weightType` | `'uniform' \| '1/cobs2' \| '1/cpred2'` | — | — | Weighting scheme |
| `modelType` | `string` | — | — | 2- or 3-compartment |
| `inputMode` | `string` | — | — | Administration route |
| `dosingParams` | `DosingParams` | — | — | Dose, tau, nDoses, tInf |
| `tEnd` | `number` | h | — | Simulation end time |

**Task output data:**

| Parameter | Type | Description |
|---|---|---|
| `bestParams` | `Map<string, number>` | Best-fit parameter values |
| `bestWssr` | `number` | WSSR at best-fit |
| `bestR2` | `number` | R² at best-fit |
| `bestAic` | `number` | AIC at best-fit |
| `bestConc` | `{t: number[], c: number[]}` | Best-fit concentration curve |
| `landscape` | `LandscapeData` | Full grid of WSSR values for visualization |
| `residuals` | `{t: number[], res: number[]}` | Residuals at best-fit |

**Computation implementation:**

| Step | Implementation method | Details | Documentation |
|---|---|---|---|
| 1 | Custom method | Generate grid: `linspace` for each selected parameter | — |
| 2 | External library in workers | `diff-grok` v1.2.0+, `mrt(task)` — solve ODE for each grid point | [diff-grok](https://github.com/datagrok-ai/diff-grok) |
| 3 | Custom method in workers | Interpolate solution at observation times, compute WSSR | — |
| 4 | Custom method | Aggregate worker results, find minimum WSSR, compute R², AIC | — |
| 5 | External library | `diff-grok` `mrt()` — solve ODE at best-fit params for full curve | — |

**Parallelization strategy:**

```
Number of workers = Math.max(1, navigator.hardwareConcurrency - 2)

Total grid points = gridResolution ^ numSelectedParams
  → distribute grid points across worker pool (round-robin or chunk)
  → each worker receives: { gridPoints[], fixedParams, modelType, inputMode, dosingParams, tEnd, observedData, weightType }
  → each worker returns: { results: { paramValues, wssr }[] }
  → coordinator aggregates, finds minimum WSSR → best-fit
```

---

### 1.2. Ports

#### Task ports: `task_primary`

| Port | Type | Interface / format | Description |
|---|---|---|---|
| Input | `PkParams` | `{ modelType, inputMode, cl, vc, vp, q, vd?, q2?, ka?, f?, dose, repeatedDosing, tau, nDoses, tInf, tEnd, mec, mtc }` | All model, dosing, and display parameters |
| Output | `PkSolution` | `{ t, ac, ap, ad?, agut?, conc, cMax, tMax, auc, cssMax?, cssMin?, t90ss? }` | Solution arrays + derived metrics |

#### Task ports: `task_generate_data`

| Port | Type | Interface / format | Description |
|---|---|---|---|
| Input | `DataGenParams` | `{ primarySolution, cv, schedule, tEnd }` | Primary solution + noise settings |
| Output | `ObservedData` | `{ tObs[], cObs[], cTrue[] }` | Sampling times and noisy observations |

#### Task ports: `task_optimize`

| Port | Type | Interface / format | Description |
|---|---|---|---|
| Input | `OptimizeParams` | `{ observedData, selectedParams[], fixedParams, gridResolution, weightType, modelType, inputMode, dosingParams, tEnd }` | Full optimization configuration |
| Output | `OptimizeResult` | `{ bestParams, bestWssr, bestR2, bestAic, bestConc, landscape, residuals }` | Best-fit results + landscape |

#### Application-level ports

| Port | Used | Description |
|---|---|---|
| Progress | Yes | `DG.TaskBarProgressIndicator` for `task_optimize` (determinate, cancelable) — reports % of grid evaluated |
| Cancellation | Yes | `pi.onCanceled` → terminates workers, resolves promises |
| Data | No | No external data loading |

### 1.3. Adapters

| Adapter | Used | Implementation |
|---|---|---|
| UI adapter | Yes | Datagrok inputs (`ui.input.float`, `ui.input.int`, `ui.input.choice`, `ui.input.bool`), `ui.bigButton`, `ui.button`, `ui.iconFA`, custom `comp_compartment_diagram`, custom `comp_opt_param_selector` |
| Display adapter | Yes | Datagrok viewers (line chart, scatter plot), custom `comp_compartment_diagram`, custom `comp_opt_landscape`, custom `comp_results_panel` |
| Worker adapter | Yes | Web Worker wrapper for `task_optimize` (`workers/optimize-worker.ts`) |
| Progress adapter | Yes | `DG.TaskBarProgressIndicator` (determinate, cancelable) for optimization |
| Data adapter | No | No external data loading |

### 1.4. Coordinator

The coordinator is the `compartmentalPkApp()` function in `app.ts`:

- **Input listening:** each input calls `debouncedRunPrimary()` from `onValueChanged`. Secondary task buttons have dedicated click handlers.
- **Reactivity management:** model type changes show/hide 3-compartment controls (V_d, Q₂); input mode changes show/hide oral controls (k_a, F) and infusion control (T_inf); repeated dosing toggle shows/hides τ and n_doses controls. See Section 9.
- **Validation trigger:** `runPrimary()` calls `validate(inputs)` before computation; validators are also attached to relevant inputs via `addValidator()` for inline hints.
- **Control state during computations:** `computationsBlocked` flag prevents `runPrimary()` during batch updates (preset loading, optimization result write, reset). Optimization button disabled during optimization.
- **Computation blocking:** `computationsBlocked = true` before batch writes, `= false` after, then single `runPrimary()`.
- **Results to display:** `updateDataFrame(result)` creates new `DG.DataFrame` and assigns to view and viewers; `updateCompartmentDiagram(result)` updates SVG; `updateMetricsPanel(result)` updates labels.
- **Resource lifecycle:** `subs[]` collects subscriptions; `activeWorkers[]` tracks workers; both cleaned up in `onViewRemoved` handler.

### 1.5. Independence Principle

Confirmed. Input controls and their reactivity (show/hide based on model type and input mode, validation hints, debounce) function independently of the computational core. The core receives a ready `PkParams` object with all parameters needed for computation; it has no knowledge of UI elements, control IDs, or conditional visibility. The `validate()` function is a pure function in `core.ts`.

---

## 2. Main View

| Field | Value |
|---|---|
| View type | `DG.TableView` |
| Description | Table view with docked concentration-time chart, compartment amounts chart, compartment diagram, and side panels for model controls (left) and optimization (right) |

---

## 3. Controls (Inputs)

### 3.1. Primary Pipeline Controls

#### Model Settings

| ID | Label | Control type | Data type | Default | Options | Nullable | Tooltip text | Group |
|---|---|---|---|---|---|---|---|---|
| `ctrl_model_type` | Model compartments | `ui.input.choice` | `string` | `'2-Compartment'` | `['2-Compartment', '3-Compartment']` | No | "Number of compartments in the PK model. Two-compartment: central + peripheral. Three-compartment adds a deep (slowly equilibrating) tissue compartment." | Model Settings |
| `ctrl_input_mode` | Administration route | `ui.input.choice` | `string` | `'IV Bolus'` | `['IV Bolus', 'IV Infusion', 'Oral']` | No | "Route of drug administration. IV Bolus: instantaneous injection. IV Infusion: constant-rate infusion over T_inf hours. Oral: first-order absorption from GI tract." | Model Settings |
| `ctrl_y_scale` | Y-axis scale | `ui.input.choice` | `string` | `'Linear'` | `['Linear', 'Semi-log']` | No | "On a semi-log scale, exponential processes appear as straight lines. Two distinct slopes = two phases: fast (α, distribution) and slow (β, elimination)." | Model Settings |

#### PK Parameters

| ID | Label | Control type | Data type | Default | Min | Max | Format | Nullable | Tooltip text | Group | Visibility condition |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `ctrl_cl` | CL — clearance (L/h) | `ui.input.float` | `number` | `5` | `0.1` | `50` | `0.0` | No | "Volume of plasma completely cleared of drug per unit time. Determines the slope of the terminal phase on a semi-log plot. Half-life t₁/₂ = 0.693 · V_ss / CL." | PK Parameters | Always |
| `ctrl_vc` | V_c — central volume (L) | `ui.input.float` | `number` | `20` | `1` | `100` | `0.0` | No | "Apparent volume of distribution in blood/plasma. Determines the initial concentration after an IV bolus: C₀ = Dose / V_c." | PK Parameters | Always |
| `ctrl_vp` | V_p — peripheral volume (L) | `ui.input.float` | `number` | `40` | `1` | `200` | `0.0` | No | "Volume of the rapidly equilibrating peripheral compartment. Larger V_p → more drug distributed to tissues → lower plasma concentrations." | PK Parameters | Always |
| `ctrl_q` | Q — intercompartmental clearance (L/h) | `ui.input.float` | `number` | `10` | `0.1` | `50` | `0.0` | No | "Rate of drug exchange between central and peripheral compartments. Large Q → fast distribution → short α-phase." | PK Parameters | Always |
| `ctrl_vd` | V_d — deep compartment volume (L) | `ui.input.float` | `number` | `60` | `1` | `200` | `0.0` | No | "Volume of the deep (slowly equilibrating) tissue compartment. Found in drugs that bind to bone, fat, or deep tissues." | PK Parameters | `ctrl_model_type = '3-Compartment'` |
| `ctrl_q2` | Q₂ — deep intercompartmental clearance (L/h) | `ui.input.float` | `number` | `2` | `0.01` | `5` | `0.00` | No | "Rate of drug exchange with the deep compartment. Typically much smaller than Q, producing a slow terminal phase." | PK Parameters | `ctrl_model_type = '3-Compartment'` |

#### Absorption Parameters

| ID | Label | Control type | Data type | Default | Min | Max | Format | Nullable | Tooltip text | Group | Visibility condition |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `ctrl_ka` | k_a — absorption rate (h⁻¹) | `ui.input.float` | `number` | `1.0` | `0.1` | `5` | `0.0` | No | "Rate of absorption from the GI tract. When k_a >> k₁₀ (no flip-flop kinetics), the peak comes quickly. When k_a ≈ k₁₀, the peak is blunted and delayed." | Absorption | `ctrl_input_mode = 'Oral'` |
| `ctrl_f` | F — bioavailability | `ui.input.float` | `number` | `0.8` | `0.01` | `1.0` | `0.00` | No | "Fraction of the dose reaching systemic circulation after oral administration. F = 1.0 for IV; typically 0.1–0.9 for oral due to first-pass metabolism." | Absorption | `ctrl_input_mode = 'Oral'` |

#### Dosing

| ID | Label | Control type | Data type | Default | Min | Max | Format | Nullable | Tooltip text | Group | Visibility condition |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `ctrl_dose` | Dose (mg) | `ui.input.float` | `number` | `100` | `1` | `1000` | `0.0` | No | "Amount of drug per dose. For IV bolus, delivered instantaneously. For infusion, delivered over T_inf hours." | Dosing | Always |
| `ctrl_repeated` | Repeated dosing | `ui.input.bool` | `boolean` | `false` | — | — | — | No | "Enable multiple-dose regimen. Shows accumulation toward steady state over repeated dosing intervals." | Dosing | Always |
| `ctrl_tau` | τ — dosing interval (h) | `ui.input.float` | `number` | `12` | `1` | `48` | `0.0` | No | "Time between successive doses. Shorter intervals → higher trough levels → faster approach to steady state." | Dosing | `ctrl_repeated = true` |
| `ctrl_n_doses` | Number of doses | `ui.input.int` | `number` | `5` | `2` | `20` | — | No | "Total number of doses administered. More doses → closer to true steady state (reached after ~5 half-lives)." | Dosing | `ctrl_repeated = true` |
| `ctrl_t_inf` | Infusion duration (h) | `ui.input.float` | `number` | `1` | `0.5` | `24` | `0.0` | No | "Duration of each constant-rate infusion. Longer infusion → lower C_max, smoother concentration profile." | Dosing | `ctrl_input_mode = 'IV Infusion'` |

#### Simulation & Display

| ID | Label | Control type | Data type | Default | Min | Max | Format | Nullable | Tooltip text | Group |
|---|---|---|---|---|---|---|---|---|---|---|
| `ctrl_t_end` | Simulation time (h) | `ui.input.float` | `number` | `72` | `1` | `240` | `0.0` | No | "Total simulation duration. For multiple dosing, ensure this covers at least n_doses × τ plus one extra interval." | Simulation |
| `ctrl_mec` | MEC — minimum effective concentration (mg/L) | `ui.input.float` | `number` | `1.0` | `0.01` | `50` | `0.00` | No | "Below this concentration, the drug has no therapeutic effect. Shown as a horizontal dashed line on the plot." | Therapeutic Window |
| `ctrl_mtc` | MTC — minimum toxic concentration (mg/L) | `ui.input.float` | `number` | `20.0` | `0.1` | `100` | `0.0` | No | "Above this concentration, toxic effects may occur. The goal of dosing is to keep concentration between MEC and MTC." | Therapeutic Window |

### 3.2. Secondary Task Triggers

| ID | Label / icon | Launches task | Tooltip text | Availability condition |
|---|---|---|---|---|
| `btn_generate_data` | `ui.bigButton('Generate Noisy Observations (CV = 15%)')` | `task_generate_data` | "Runs the model at current parameters, samples at clinical time points, and adds log-normal noise. Try recovering the true parameters with the optimizer — higher CV makes it harder!" | Always (label updates dynamically with current CV value) |
| `btn_run_optimize` | `ui.bigButton('Run Grid Search (estimated: 400 ODE solves)')` | `task_optimize` | "Perform exhaustive search over selected parameters. Each grid point requires a full ODE solve. The parameter combination minimizing WSSR is returned as the best fit." | Requires observed data from `task_generate_data`; disabled during optimization |

### 3.3. Secondary Task Controls

#### Task controls: `task_generate_data`

UI type: inline controls in optimization panel (not a dialog).

| ID | Label | Control type | Data type | Default | Min | Max | Format | Nullable | Tooltip text |
|---|---|---|---|---|---|---|---|---|---|
| `ctrl_noise_cv` | Observation noise CV (%) | `ui.input.float` | `number` | `15` | `5` | `40` | `0.0` | No | "Coefficient of variation for log-normal noise. Higher CV = noisier data = harder optimization. 15% is typical for clinical PK data." |
| `ctrl_sampling` | Sampling schedule | `ui.input.choice` | `string` | `'Standard PK'` | — | — | — | No | "Time points at which samples are taken. Standard PK: 7 points. Rich: 15 points (better identifiability). Sparse: 4 points (realistic clinical limitation)." |

Options for `ctrl_sampling`: `['Standard PK', 'Rich', 'Sparse']`

#### Task controls: `task_optimize`

UI type: inline controls in optimization panel.

| ID | Label | Control type | Data type | Default | Min | Max | Format | Nullable | Tooltip text |
|---|---|---|---|---|---|---|---|---|---|
| `ctrl_grid_res` | Grid density per axis | `ui.input.int` | `number` | `20` | `5` | `50` | — | No | "Number of points per parameter axis. Total ODE solves = grid_points ^ number_of_parameters. For 2 parameters at 30 grid points = 900 solves — usually takes < 2 seconds with MRT." |
| `ctrl_weight` | Objective function weighting | `ui.input.choice` | `string` | `'1/C_obs²'` | — | — | — | No | "WSSR weighting scheme. 1/C_obs² is standard in PK because measurement error is typically proportional to concentration. Uniform weights treat all residuals equally." |

Options for `ctrl_weight`: `['Uniform', '1/C_obs²', '1/C_pred²']`

### 3.4. Other Buttons and Actions

| ID | Label / icon | Action | Tooltip text | Availability condition |
|---|---|---|---|---|
| `btn_preset_vanc` | `ui.button('Vancomycin (CL=4.6, V_c=28, t₁/₂≈6h)')` | Load vancomycin PK parameters: CL=4.6, V_c=28, V_p=28, Q=8.8, Dose=1000, 2-compartment, IV Infusion, T_inf=1 | "Load real-world PK parameters for vancomycin, a glycopeptide antibiotic with two-compartment kinetics." | Always |
| `btn_preset_theo` | `ui.button('Theophylline (CL=3.5, V_c=30, t₁/₂≈8h)')` | Load theophylline PK parameters: CL=3.5, V_c=30, V_p=20, Q=5, Dose=300, 2-compartment, Oral, k_a=1.5, F=0.9 | "Load real-world PK parameters for theophylline, a bronchodilator with oral administration and well-characterized PK." | Always |
| `btn_preset_metf` | `ui.button('Metformin (CL=26, V_c=63, t₁/₂≈5h)')` | Load metformin PK parameters: CL=26, V_c=63, V_p=100, Q=15, Dose=500, 2-compartment, Oral, k_a=2.0, F=0.55 | "Load real-world PK parameters for metformin, a diabetes drug with high clearance and moderate bioavailability." | Always |
| `btn_reset` | `ui.iconFA('undo')` | Reset all controls to default values | "Reset all parameters to default values" | Always |
| `btn_compare` | Compare: sliders vs optimized | `ui.input.bool` | Toggle overlay of current slider curve and best-fit curve on concentration-time chart | "Overlay both curves: the current slider-set parameters and the optimized best-fit. See how they differ." | Requires optimization results |

### 3.5. Custom UI Components

| Component ID | Brief description | Role | UI component specification |
|---|---|---|---|
| `comp_compartment_diagram` | Animated SVG compartment diagram. Rectangles labeled "Central", "Peripheral", "Deep" (if 3-comp) with directional arrows. Arrow thickness proportional to current mass flow at time t. Arrow labels show rate constants with physiological values (e.g., "CL = 5 L/h"). Includes time slider and play/pause button for animation. Adapts layout for 2-comp vs 3-comp. | Display | See description below |
| `comp_opt_param_selector` | Parameter selection panel for optimization. Shows checkboxes for available PK parameters (CL, V_c, V_p, Q, and conditionally V_d, Q₂, k_a, F). For each checked parameter, shows min/max range sliders (default ±50% of current value). Shows estimated total ODE solves count. Enforces 1–3 parameter limit. | Control | See description below |
| `comp_results_panel` | Display panel showing best-fit parameters, goodness-of-fit metrics (WSSR, R², AIC). Built with `ui.divV`/`ui.label`; values updated via `textContent`. | Display | Inline |
| `comp_pk_metrics` | Display panel showing PK metrics: C_max, t_max, AUC, and steady-state metrics when applicable. Built with `ui.divV`/`ui.label`. | Display | Inline |

**`comp_compartment_diagram` description:**

- **Visual elements:** Rounded rectangles for compartments (Central, Peripheral, Deep). Arrows between compartments with arrowheads. Arrow thickness ∝ mass flow = rate_constant × amount_in_source. Labels on arrows show both the rate constant symbol and current value. Elimination arrow from Central pointing down/out.
- **States:** 2-compartment (2 boxes, 3 arrows: k₁₂, k₂₁, k₁₀) vs 3-compartment (3 boxes, 5 arrows: k₁₂, k₂₁, k₁₃, k₃₁, k₁₀). For oral: additional "Gut" box with k_a arrow to Central.
- **Animation controls:** Time slider (0 to T_end), play/pause button, speed control. When playing, time advances and arrow thicknesses update.
- **Styles:** CSS prefix `cpk-diagram-`.

**`comp_opt_param_selector` description:**

- **Visual elements:** Checkbox + label for each parameter. Below each checked parameter: min/max range inputs with current value shown. Summary line: "Total ODE solves: {N}" with warning icon if N > 10,000.
- **States:** Parameters dynamically available based on model type and input mode. Checked parameters limited to 1–3 (further checkboxes disabled when 3 are selected).
- **Events:** Emits `onSelectionChanged` with `{ selectedParams: {name, min, max}[] }`.
- **Styles:** CSS prefix `cpk-opt-selector-`.

---

## 4. Result Display Elements

### 4.1. Primary Pipeline Display Elements

| ID | Type | Associated output data | Docking location |
|---|---|---|---|
| `view_conc_time` | Datagrok viewer `line chart` | `Time`, `Concentration` (+ MEC/MTC lines + therapeutic window shading on linear scale) | Main area, top-right |
| `view_compartments` | Datagrok viewer `line chart` | `Time`, `A_central`, `A_peripheral`, [`A_deep`], [`A_gut`] | Main area, bottom-left |
| `comp_compartment_diagram` | Custom HTMLElement (SVG) | `t`, `ac`, `ap`, [`ad`], [`agut`] — animated | Main area, bottom-right |
| `comp_pk_metrics` | Custom HTMLElement | `cMax`, `tMax`, `auc`, `cssMax`, `cssMin`, `t90ss` | Left panel, below controls |
| `view_table` | `DG.TableView` grid | All columns | Center (default) |

**Details for `view_conc_time`:**

| Property | Value |
|---|---|
| X axis | `Time (h)` |
| Y axis | `Concentration (mg/L)` — linear or log scale based on `ctrl_y_scale` |
| Series | Concentration — solid line |
| Markers | C_max point (labeled), C_ss_max point, C_ss_min point, t_90_ss point (when applicable) |
| Horizontal lines | MEC (dashed, labeled), MTC (dashed, labeled) |
| Shading | Therapeutic window (area between MEC and MTC) on linear scale |
| Overlay (when observed data exists) | Scatter points for C_obs |
| Overlay (when optimization done + compare on) | Best-fit curve as distinct line style |

**Details for `view_compartments`:**

| Property | Value |
|---|---|
| X axis | `Time (h)` |
| Y axis | `Amount (mg)` |
| Series | A_central (solid), A_peripheral (solid), A_deep (solid, 3-comp only), A_gut (solid, oral only) |
| Legend | Full compartment names |

### 4.2. Secondary Task Display Elements

#### Task display: `task_generate_data`

Observed data points are overlaid on `view_conc_time` as scatter points (circles). No separate display element — integrated into the primary line chart.

| ID | Type | Associated output data | Placement |
|---|---|---|---|
| (integrated into `view_conc_time`) | Scatter overlay | `tObs`, `cObs` | Overlaid on concentration-time chart |

#### Task display: `task_optimize`

| ID | Type | Associated output data | Placement |
|---|---|---|---|
| `comp_results_panel` | Custom HTMLElement | `bestParams`, `bestWssr`, `bestR2`, `bestAic` | Right panel, optimization section |
| `view_landscape` | Datagrok viewer or custom | `landscape` data | Right panel, below results |
| `view_residuals` | Datagrok viewer `scatter plot` | `residuals.t`, `residuals.res` | Right panel, below landscape |

**Details for `view_landscape`:**

| Number of optimized params | Visualization |
|---|---|
| 1 | Line chart: WSSR vs parameter value. Minimum marked with vertical line. |
| 2 | Heatmap: WSSR over 2D parameter grid. Minimum marked with crosshair. Color: lower WSSR = better. |
| 3 | Three 2D heatmaps (pairwise projections). Each shows minimum WSSR marginalized over the third parameter. |

**Details for `view_residuals`:**

| Property | Value |
|---|---|
| X axis | `Time (h)` |
| Y axis | `C_obs − C_pred (mg/L)` |
| Horizontal line | y = 0 (reference) |
| Series | Residuals for current slider params (if compare on) + residuals for best-fit |
| Title | "Residuals — Observed minus Predicted" |

---

## 5. Layout and UI Element Placement

### 5.1. Control Placement

| Area | Content (control IDs) |
|---|---|
| Left panel (top section) | Model Settings: `ctrl_model_type`, `ctrl_input_mode`, `ctrl_y_scale` |
| Left panel (PK section) | PK Parameters: `ctrl_cl`, `ctrl_vc`, `ctrl_vp`, `ctrl_q`, `ctrl_vd`, `ctrl_q2` |
| Left panel (absorption section) | Absorption: `ctrl_ka`, `ctrl_f` |
| Left panel (dosing section) | Dosing: `ctrl_dose`, `ctrl_repeated`, `ctrl_tau`, `ctrl_n_doses`, `ctrl_t_inf` |
| Left panel (sim section) | Simulation: `ctrl_t_end` |
| Left panel (window section) | Therapeutic Window: `ctrl_mec`, `ctrl_mtc` |
| Left panel (metrics section) | `comp_pk_metrics` |
| Ribbon (Presets group) | `btn_preset_vanc`, `btn_preset_theo`, `btn_preset_metf` |
| Ribbon (Actions group) | `btn_reset` |
| Right panel (data gen section) | `ctrl_noise_cv`, `ctrl_sampling`, `btn_generate_data` |
| Right panel (opt setup section) | `comp_opt_param_selector`, `ctrl_grid_res`, `ctrl_weight`, `btn_run_optimize` |
| Right panel (results section) | `comp_results_panel`, `btn_compare` |
| Right panel (landscape section) | `view_landscape` |
| Right panel (residuals section) | `view_residuals` |

**Structure of left panel form:**

```
ui.h2('Model Settings')
  ctrl_model_type
  ctrl_input_mode
  ctrl_y_scale

ui.h2('PK Parameters')
  ctrl_cl
  ctrl_vc
  ctrl_vp
  ctrl_q
  ctrl_vd  (shown if 3-Compartment)
  ctrl_q2  (shown if 3-Compartment)

ui.h2('Absorption')  (shown if Oral)
  ctrl_ka
  ctrl_f

ui.h2('Dosing')
  ctrl_dose
  ctrl_repeated
  ctrl_tau  (shown if repeated)
  ctrl_n_doses  (shown if repeated)
  ctrl_t_inf  (shown if IV Infusion)

ui.h2('Simulation')
  ctrl_t_end

ui.h2('Therapeutic Window')
  ctrl_mec
  ctrl_mtc

ui.h2('PK Metrics')
  comp_pk_metrics
```

**Structure of right panel (Brute-Force Parameter Estimation):**

```
ui.h2('Synthetic Observation Data')
  ctrl_noise_cv
  ctrl_sampling
  btn_generate_data

ui.h2('Optimization Setup')
  comp_opt_param_selector
  ctrl_grid_res
  ctrl_weight
  btn_run_optimize

ui.h2('Best-Fit Parameters')  (shown after optimization)
  comp_results_panel
  btn_compare

view_landscape  (shown after optimization)
view_residuals  (shown after optimization)
```

### 5.2. Display Element Placement

| Element ID | Docking area | Position / ratio |
|---|---|---|
| Left panel form | `DG.DOCK_TYPE.LEFT` relative to root | ratio `0.3` |
| Right panel (optimization) | `DG.DOCK_TYPE.RIGHT` relative to root | ratio `0.3` |
| `view_conc_time` | `DG.DOCK_TYPE.TOP` relative to center area | ratio `0.5` |
| `view_compartments` | `DG.DOCK_TYPE.LEFT` relative to bottom center | ratio `0.5` |
| `comp_compartment_diagram` | `DG.DOCK_TYPE.RIGHT` relative to bottom center | ratio `0.5` |
| `view_table` (grid) | Default position (center) | — |

### 5.3. Styles

CSS file: `css/compartmental-pk.css`. Import: `import '../../css/compartmental-pk.css'`. All custom classes use the `cpk-` prefix for isolation.

**Static styles:**

| Element | CSS class(es) | Description |
|---|---|---|
| PK metrics panel | `.cpk-metrics-panel` | Background, padding, border-radius, gap |
| Metrics value | `.cpk-metrics-value` | Monospace font, right-aligned |
| Results panel | `.cpk-results-panel` | Background, padding for best-fit display |
| Results value | `.cpk-results-value` | Monospace font, right-aligned |
| Diagram container | `.cpk-diagram-container` | SVG container with border |
| Diagram compartment box | `.cpk-diagram-compartment` | Rounded rect, fill, stroke |
| Diagram arrow | `.cpk-diagram-arrow` | Stroke, marker-end |
| Diagram label | `.cpk-diagram-label` | Font size, alignment |
| Optimization section headers | `.cpk-opt-section-header` | Bold, separator |
| Parameter selector | `.cpk-opt-selector-row` | Flex row, checkbox + slider layout |
| Warning text | `.cpk-warning` | Orange color, italic |

**Dynamic styles:**

| Element | CSS class(es) | Condition | Description |
|---|---|---|---|
| Optimize button | `.cpk-btn--disabled` | During optimization or no observed data | Disables pointer events, reduces opacity |
| Generate data button | `.cpk-btn--active` | Data has been generated | Highlighted border to indicate data exists |
| Conditional controls | `.cpk-hidden` | When control is not applicable | `display: none` for conditional show/hide |
| Diagram arrow | `stroke-width` (inline) | Proportional to mass flow | Dynamic SVG attribute, not CSS class |

---

## 6. User Feedback

### 6.1. Control Tooltips

All tooltips are defined in Section 3 (Tooltip text column). Mechanism: `tooltipText` property for `ui.input.*` controls.

| Control ID | Tooltip text | Mechanism |
|---|---|---|
| `ctrl_model_type` | (see Section 3.1) | `tooltipText` property |
| `ctrl_input_mode` | (see Section 3.1) | `tooltipText` property |
| `ctrl_y_scale` | (see Section 3.1) | `tooltipText` property |
| `ctrl_cl` | (see Section 3.1) | `tooltipText` property |
| `ctrl_vc` | (see Section 3.1) | `tooltipText` property |
| `ctrl_vp` | (see Section 3.1) | `tooltipText` property |
| `ctrl_q` | (see Section 3.1) | `tooltipText` property |
| `ctrl_vd` | (see Section 3.1) | `tooltipText` property |
| `ctrl_q2` | (see Section 3.1) | `tooltipText` property |
| `ctrl_ka` | (see Section 3.1) | `tooltipText` property |
| `ctrl_f` | (see Section 3.1) | `tooltipText` property |
| `ctrl_dose` | (see Section 3.1) | `tooltipText` property |
| `ctrl_repeated` | (see Section 3.1) | `tooltipText` property |
| `ctrl_tau` | (see Section 3.1) | `tooltipText` property |
| `ctrl_n_doses` | (see Section 3.1) | `tooltipText` property |
| `ctrl_t_inf` | (see Section 3.1) | `tooltipText` property |
| `ctrl_t_end` | (see Section 3.1) | `tooltipText` property |
| `ctrl_mec` | (see Section 3.1) | `tooltipText` property |
| `ctrl_mtc` | (see Section 3.1) | `tooltipText` property |
| `ctrl_noise_cv` | (see Section 3.3) | `tooltipText` property |
| `ctrl_sampling` | (see Section 3.3) | `tooltipText` property |
| `ctrl_grid_res` | (see Section 3.3) | `tooltipText` property |
| `ctrl_weight` | (see Section 3.3) | `tooltipText` property |
| `btn_generate_data` | (see Section 3.2) | `ui.bigButton` third argument |
| `btn_run_optimize` | (see Section 3.2) | `ui.bigButton` third argument |
| `btn_reset` | "Reset all parameters to default values" | `ui.iconFA` third argument |
| `btn_preset_vanc` | (see Section 3.4) | `ui.button` tooltip |
| `btn_preset_theo` | (see Section 3.4) | `ui.button` tooltip |
| `btn_preset_metf` | (see Section 3.4) | `ui.button` tooltip |
| `btn_compare` | (see Section 3.4) | `tooltipText` property |

### 6.2. Validators as Feedback

| Input ID | Validation source | Description |
|---|---|---|
| `ctrl_cl` … `ctrl_t_end`, `ctrl_mec`, `ctrl_mtc` | `validate()` from core | Each input has `addValidator()` calling `validate(getInputs())` and returning the error for its own ID, or `null`. Displays inline hint on invalid value. |

### 6.3. Progress Bar

| Task | Progress bar | Type | Cancellation support |
|---|---|---|---|
| `task_primary` | No | — | — |
| `task_generate_data` | No | — | — |
| `task_optimize` | Yes | Determinate (`DG.TaskBarProgressIndicator`) — reports % of grid points evaluated | Yes (`pi.onCanceled` → terminates workers) |

---

## 7. Validation

### 7.1. Primary Pipeline Validation

#### Complex Validation Rules

| Rule ID | Condition (invalid) | Affected inputs (ID) | Error message |
|---|---|---|---|
| `val_01` | `cl ≤ 0` | `ctrl_cl` | "Clearance must be positive" |
| `val_02` | `vc ≤ 0` | `ctrl_vc` | "Central volume must be positive" |
| `val_03` | `vp ≤ 0` | `ctrl_vp` | "Peripheral volume must be positive" |
| `val_04` | `q ≤ 0` | `ctrl_q` | "Intercompartmental clearance must be positive" |
| `val_05` | `vd ≤ 0` (when 3-comp) | `ctrl_vd` | "Deep compartment volume must be positive" |
| `val_06` | `q2 ≤ 0` (when 3-comp) | `ctrl_q2` | "Deep intercompartmental clearance must be positive" |
| `val_07` | `ka ≤ 0` (when oral) | `ctrl_ka` | "Absorption rate must be positive" |
| `val_08` | `f ≤ 0 or f > 1` (when oral) | `ctrl_f` | "Bioavailability must be between 0 (exclusive) and 1 (inclusive)" |
| `val_09` | `dose ≤ 0` | `ctrl_dose` | "Dose must be positive" |
| `val_10` | `tau ≤ 0` (when repeated) | `ctrl_tau` | "Dosing interval must be positive" |
| `val_11` | `nDoses < 1` (when repeated) | `ctrl_n_doses` | "Number of doses must be at least 1" |
| `val_12` | `tInf ≤ 0` (when infusion) | `ctrl_t_inf` | "Infusion duration must be positive" |
| `val_13` | `tInf > tau` (when infusion + repeated) | `ctrl_t_inf` | "Infusion duration cannot exceed dosing interval" |
| `val_14` | `tEnd ≤ 0` | `ctrl_t_end` | "Simulation time must be positive" |
| `val_15` | `mtc ≤ mec` | `ctrl_mtc` | "MTC must be greater than MEC" |
| `val_16` | `mec < 0` | `ctrl_mec` | "MEC cannot be negative" |

#### Validation Order

```
1. val_01 through val_12 (independent, all checked)
2. val_13 (checked only if val_10 and val_12 both passed — needs valid tau and tInf)
3. val_14 through val_16 (independent)
```

#### Returned Map Format

```
Map<InputId, string>
```

Where `InputId` is the control ID string (e.g., `'ctrl_cl'`).

### 7.2. Secondary Task Validation

#### Task validation: `task_optimize`

| Rule ID | Condition (invalid) | Affected inputs (ID) | Error message |
|---|---|---|---|
| `val_opt_01` | No observed data available | — | "Generate synthetic observations first" |
| `val_opt_02` | No parameters selected for optimization | — | "Select at least one parameter to optimize" |
| `val_opt_03` | More than 3 parameters selected | — | "Select at most 3 parameters for optimization" |
| `val_opt_04` | Any selected param has min ≥ max | `comp_opt_param_selector` | "Search range minimum must be less than maximum for {param}" |
| `val_opt_05` | Primary validation fails | — | "Fix primary parameter errors before optimizing" |

Validation order:
```
1. val_opt_05 (primary params must be valid)
2. val_opt_01 (observed data must exist)
3. val_opt_02, val_opt_03, val_opt_04 (optimization config)
```

Return format: `{ errors: Map<string, string>, warning: string | null }`

Warning when `gridResolution^numParams > 10000`: "This may take a few seconds. Total ODE solves: {N}"

---

## 8. Main Pipeline

### 8.1. Primary Pipeline

| Step | Description |
|---|---|
| Parameter input | User sets values through controls → `getInputs()` converts to `PkParams` (respecting current model type and input mode) |
| Validation | `validate(inputs)` returns `Map<InputId, string>`. On errors: mark invalid inputs, keep previous results or clear |
| Computation | `solvePk(inputs)` — synchronous ODE solution via segmented `mrt()` calls |
| Result display | `updateDataFrame(result)` creates new DataFrame; `updateChart()` updates line chart; `updateCompartmentDiagram(result)` updates SVG; `updateMetrics(result)` updates labels |

Reactive trigger: `onValueChanged` with debounce 100 ms.

Error behavior on validation failure: keep previous results displayed, show validation hints on affected inputs.

Error behavior on computation failure: clear results + `grok.shell.error(msg)`.

### 8.2. Secondary Pipelines

#### Task pipeline: `task_generate_data`

| Step | Description |
|---|---|
| Trigger | `btn_generate_data` click |
| Custom UI | N/A — controls are inline in optimization panel |
| Validation | Primary validation must pass; noise CV and sampling schedule have their own constraints |
| Computation | Interpolate primary solution at sampling times, add log-normal noise |
| Result display | Observed data points shown as scatter overlay on `view_conc_time` |
| Feedback to primary | No — observed data is stored separately, does not affect primary controls |

#### Task pipeline: `task_optimize`

| Step | Description |
|---|---|
| Trigger | `btn_run_optimize` click |
| Custom UI | N/A — controls are inline in optimization panel |
| Validation | `validateOptimization()` checks: observed data exists, 1–3 params selected, valid ranges, primary params valid |
| Computation | Generate grid → distribute across WebWorkers → collect results → find min WSSR → compute best-fit curve |
| Result display | `comp_results_panel` shows best-fit params and metrics; `view_landscape` shows WSSR landscape; `view_residuals` shows residual plot |
| Feedback to primary | When `btn_compare` is toggled on: best-fit curve overlaid on `view_conc_time`. Optionally, user can click "Apply Best-Fit" to write best-fit params to primary controls (batch update). |

### 8.3. Common Pipeline Aspects

#### Control Behavior During Computations

| Pipeline / task | Controls blocked | Which controls |
|---|---|---|
| `task_primary` | No | — (computation < 100 ms) |
| `task_generate_data` | No | — (computation < 10 ms) |
| `task_optimize` | `btn_run_optimize` only | `btn_run_optimize` disabled via `.cpk-btn--disabled`; all other controls remain interactive (user can adjust sliders while optimization runs) |

#### Computation Error Handling

| Pipeline / task | Strategy | Notification method |
|---|---|---|
| `task_primary` | Clear results | `grok.shell.error` |
| `task_generate_data` | No data generated, show message | `grok.shell.error` |
| `task_optimize` (partial failures) | Skip failed grid points, use valid results | `grok.shell.warning` |
| `task_optimize` (all fail) | Abort optimization | `grok.shell.error` |

### 8.4. Computation Blocking and Batch Input Updates

| Scenario | Source | Target controls (ID) | Blocked pipelines |
|---|---|---|---|
| Preset loading | `btn_preset_*` | `ctrl_model_type`, `ctrl_input_mode`, `ctrl_cl`, `ctrl_vc`, `ctrl_vp`, `ctrl_q`, `ctrl_dose`, `ctrl_ka`, `ctrl_f`, `ctrl_t_inf` | Primary (blocked) |
| Apply best-fit | `task_optimize` results (optional action) | Selected optimized params (e.g., `ctrl_cl`, `ctrl_vc`) | Primary (blocked) |
| Reset to defaults | `btn_reset` | All `ctrl_*` | Primary (blocked) |
| Format setting on init | Coordinator | All numeric `ctrl_*` | Primary (blocked) |

Reactivity mode during batch update:

| Scenario | Reactivity mode |
|---|---|
| All above | `computationsBlocked = true` → writes → `computationsBlocked = false` → single `runPrimary()` |

---

## 9. Reactivity and Dependencies Between Inputs

### 9.1. Dependency Graph

| Source (input ID) | Target (input IDs) | Reaction type | Logic |
|---|---|---|---|
| `ctrl_model_type` | `ctrl_vd`, `ctrl_q2` | Availability | Show if '3-Compartment', hide otherwise |
| `ctrl_model_type` | `comp_opt_param_selector` | Option list | Update available parameters for optimization |
| `ctrl_model_type` | `comp_compartment_diagram` | Display update | Switch between 2-comp and 3-comp diagram layout |
| `ctrl_input_mode` | `ctrl_ka`, `ctrl_f` | Availability | Show if 'Oral', hide otherwise |
| `ctrl_input_mode` | `ctrl_t_inf` | Availability | Show if 'IV Infusion', hide otherwise |
| `ctrl_input_mode` | `comp_opt_param_selector` | Option list | Add/remove k_a, F from available optimization parameters |
| `ctrl_input_mode` | `comp_compartment_diagram` | Display update | Show/hide Gut compartment |
| `ctrl_repeated` | `ctrl_tau`, `ctrl_n_doses` | Availability | Show if true, hide otherwise |
| `ctrl_noise_cv` | `btn_generate_data` | Label update | Update button label: "Generate Noisy Observations (CV = {cv}%)" |
| `ctrl_grid_res` | `btn_run_optimize` | Label update | Update button label with estimated ODE solves count |
| `comp_opt_param_selector` | `btn_run_optimize` | Label update | Update estimated solves: `gridRes^numParams` |
| Any primary `ctrl_*` | Primary pipeline | Reactive trigger | Rerun `task_primary` |

### 9.2. Debounce / Throttle

| Input ID | Strategy | Interval (ms) |
|---|---|---|
| All numeric `ctrl_*` (float/int) | debounce | 100 |
| `ctrl_model_type`, `ctrl_input_mode`, `ctrl_y_scale`, `ctrl_repeated`, `ctrl_sampling` | immediate | 0 (no debounce for choice/toggle) |

---

## 10. Data Lifecycle

### 10.1. Data Input

Primary method: manual input via left panel controls.

Initial state: all controls initialized with defaults. `task_primary` runs automatically on initialization.

### 10.2. Loading from Resources

Preset buttons load predefined parameter sets (see Section 3.4). No external file or network loading.

| Trigger (button ID) | Resource | Format | Mapping to inputs (ID) |
|---|---|---|---|
| `btn_preset_vanc` | Hardcoded constants | Object literal | `ctrl_cl=4.6, ctrl_vc=28, ctrl_vp=28, ctrl_q=8.8, ctrl_dose=1000, ctrl_model_type='2-Compartment', ctrl_input_mode='IV Infusion', ctrl_t_inf=1` |
| `btn_preset_theo` | Hardcoded constants | Object literal | `ctrl_cl=3.5, ctrl_vc=30, ctrl_vp=20, ctrl_q=5, ctrl_dose=300, ctrl_model_type='2-Compartment', ctrl_input_mode='Oral', ctrl_ka=1.5, ctrl_f=0.9` |
| `btn_preset_metf` | Hardcoded constants | Object literal | `ctrl_cl=26, ctrl_vc=63, ctrl_vp=100, ctrl_q=15, ctrl_dose=500, ctrl_model_type='2-Compartment', ctrl_input_mode='Oral', ctrl_ka=2.0, ctrl_f=0.55` |

### 10.3. Results Table Lifecycle

Update strategy: **DataFrame replacement** — on each recomputation a new `DG.DataFrame` is created and assigned to `view.dataFrame` and all viewers.

```
1. Application initialization
   → solvePk with DEFAULTS → new DG.DataFrame with columns [Time, Concentration, A_central, A_peripheral, (A_deep), (A_gut)]
   → DataFrame added to TableView

2. task_primary completion
   → new DataFrame created from solution arrays
   → assigned to view.dataFrame and all viewer dataFrames
   → All viewers redrawn reactively

3. task_generate_data completion
   → Observed data stored in application state (not in main DataFrame)
   → Scatter overlay added/updated on view_conc_time

4. task_optimize completion
   → Optimization results stored in application state
   → Results panel, landscape, and residuals displays updated
   → Best-fit curve stored for optional overlay

5. btn_reset press
   → All controls reset to defaults → task_primary runs → new DataFrame
   → Observed data and optimization results cleared
```

---

## 11. Error Handling Beyond Computations

| Error type | Strategy | Notification method |
|---|---|---|
| ODE solver error (numerical instability) | Clear results, show message | `grok.shell.error` |
| Worker creation error | Abort optimization, show message | `grok.shell.error` |
| Partial worker errors | Skip failed grid points, use valid results | `grok.shell.warning` with count of failures |
| All workers fail | Abort optimization | `grok.shell.error` |
| Interpolation error (sample time outside solution range) | Skip affected sample points, warn | `grok.shell.warning` |

---

## 12. Subscriptions and Resource Management

### 12.1. Event Subscriptions

| Subscription | Event | Cleanup mechanism |
|---|---|---|
| View removal | `grok.events.onViewRemoved` | `sub.unsubscribe()` in cleanup handler |
| Input change (×N inputs) | `input.onValueChanged` | Collected in `subs[]`, unsubscribed on close |
| Optimization cancel | `pi.onCanceled` | `cancelSub.unsubscribe()` after workers complete |
| Diagram animation timer | `setInterval` / `requestAnimationFrame` | `clearInterval` / `cancelAnimationFrame` on close |

All subscriptions are collected in `subs[]` and unsubscribed when the view is removed.

### 12.2. Worker Termination

| Worker pool | Created in | Termination mechanism |
|---|---|---|
| Optimization workers | `runOptimization()` | `w.terminate()` for each worker in `terminateWorkers()`, called on completion, cancellation, and view close |

---

## 13. Application Closure

On view close (`grok.events.onViewRemoved`), the coordinator performs:

- [x] All event subscriptions unsubscribed (`subs[]` iteration)
- [x] All web workers terminated (`terminateWorkers()`)
- [x] Pending debounce timer cleared (`clearTimeout(debounceTimer)`)
- [x] Diagram animation timer cleared
- [x] No open secondary task dialogs to close (optimization uses inline controls)

Closure handler: `grok.events.onViewRemoved.subscribe(...)` — checks `v.name` matches application view name.

---

## 14. Accessibility and UX

### 14.1. Keyboard Shortcuts

N/A — no custom keyboard shortcuts.

### 14.2. Context Menus

N/A — no custom context menus.

### 14.3. Undo / Redo

Not supported. The `btn_reset` button resets all controls to default values. Disease presets provide quick parameter loading as an alternative.

---

## 15. Testing

### 15.1. Computational Part (Core)

| File | Categories | Test count | Description |
|---|---|---|---|
| `src/tests/pk-math-tests.ts` | Math: PK ODE RHS, Math: PK Solve, Math: PK Metrics, Math: Data Generation, Math: Optimization | ~30 | ODE RHS verification for all model/route combinations, solve properties, PK metrics, noise model, WSSR computation |
| `src/tests/pk-validation-tests.ts` | API: Validation | ~25 | Validation rules for all inputs, conditional rules, boundary cases, multiple errors |

Tests are run via `grok test` (entry point: `src/package-test.ts`).

### 15.2. Inputs

| Category | Coverage | Description |
|---|---|---|
| Boundary values | Lower/upper bounds for all numeric controls | e.g., `cl=0.1`, `f=1.0`, `nDoses=1` |
| Invalid values | Zero and negative for positive-domain params | e.g., `cl=0`, `cl=-1`, `f=0`, `f=1.1` |
| Invalid combinations | Cross-parameter constraints | e.g., `tInf > tau`, `mtc ≤ mec` |
| Conditional validation | Rules that apply only for certain model/mode | e.g., `vd ≤ 0` only checked when 3-compartment |
| Multiple simultaneous errors | Several invalid inputs at once | e.g., `cl=0, vc=0, dose=0` → ≥ 3 errors |
| Valid defaults | All defaults pass validation | DEFAULTS → `errors.size = 0` |

### 15.3. Mathematical Verification

#### Level 1 Verification (required)

**Formula/equation verification:**

| Test category | Test count | What is verified | Reference source |
|---|---|---|---|
| Math: PK ODE RHS (2-comp, IV Bolus) | 2 | ODE RHS at t=0 and t>0 for 2-comp bolus | Manual calculation (see Ref #1 in Section 1.1) |
| Math: PK ODE RHS (3-comp, IV Bolus) | 2 | ODE RHS at t=0 for 3-comp bolus | Manual calculation (see Ref #2 in Section 1.1) |
| Math: PK ODE RHS (2-comp, Oral) | 2 | ODE RHS at t=0 for oral with gut compartment | Manual calculation (see Ref #3 in Section 1.1) |
| Math: PK ODE RHS (2-comp, IV Infusion) | 2 | ODE RHS during and after infusion | Manual calculation (see Ref #4 in Section 1.1) |
| Math: Micro-constants | 3 | k₁₀, k₁₂, k₂₁, k₁₃, k₃₁ from physiological params | Manual: k₁₀=CL/V_c, k₁₂=Q/V_c, etc. |
| Math: WSSR computation | 2 | WSSR with known inputs for each weight type | Manual calculation |

**Output property verification:**

| Test category | Test count | Properties verified |
|---|---|---|
| Math: PK Solve properties | 8 | Non-negative amounts, initial conditions (C(0)=Dose/V_c for bolus, C(0)=0 for others), mass balance, decay to zero, C_max > 0, t_max > 0, steady-state max > min |

#### Level 2 Verification

**Numerical method verification:**

| Test category | Test count | Reference problems | Tolerance | Source |
|---|---|---|---|---|
| Math: MRT solver | 2 | Bi-exponential analytical solution for 2-comp IV bolus; comparison of numerical vs analytical C(t) at multiple time points | Max relative error < 1% | Analytical solution: C(t) = A·exp(−αt) + B·exp(−βt) from eigenvalue decomposition |

**Convergence verification:**

| Status | Description |
|---|---|
| Not yet covered | Could verify that reducing MRT tolerance produces converging C(t) profiles |

**Asymptotic/equilibrium behavior:**

| Status | Description |
|---|---|
| Implemented | Verify C(T_end) → 0 for large T_end (single dose); verify C_ss convergence for multiple dosing |

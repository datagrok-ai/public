# Application Specification: SIR / SEIR Epidemic Simulation

<!-- TOC for partial reading -->
<!-- Section 1: General Architecture — lines 15–250 -->
<!-- Section 2: Main View — lines 251–259 -->
<!-- Section 3: Controls — lines 260–345 -->
<!-- Section 4: Display Elements — lines 346–420 -->
<!-- Section 5: Layout — lines 421–500 -->
<!-- Section 6: User Feedback — lines 501–545 -->
<!-- Section 7: Validation — lines 546–585 -->
<!-- Section 8: Pipeline — lines 586–650 -->
<!-- Section 9: Reactivity — lines 651–675 -->
<!-- Sections 10–15: Data, Errors, Resources, Closure, UX, Testing — lines 676–900 -->

## 1. General Architecture

### 1.0. General Information

| Field | Value |
|---|---|
| Application name | SIR / SEIR Epidemic Simulation |
| Package | SirSeirEpidemic |
| Entry function | `sirSeirEpidemicApp()` |
| Brief description | Interactive ODE simulation of SIR and SEIR epidemiological models with epidemic dynamics, effective reproduction number, phase portrait, disease presets, and vaccination coverage. |
| Main view type | `DG.TableView` |

### 1.1. Core

#### Task List

| Task ID | Name | Pipeline type | Trigger | Synchronicity | Execution environment | Parallelization |
|---|---|---|---|---|---|---|
| `task_primary` | SIR/SEIR epidemic simulation | Primary (reactive) | Any input change | Sync | Main thread | No |

No secondary tasks. Disease presets are batch updates to primary controls, not separate computations.

#### Dependencies Between Tasks

```
task_primary — independent (single task)
```

#### Task: `task_primary`

**General characteristics:**

| Field | Value |
|---|---|
| Name | SIR/SEIR epidemic simulation |
| Description | Numerical solution of the SIR or SEIR ODE system using MRT, returning time series for all compartments, R_eff(t), peak infection data, and summary statistics |
| Dependency on other tasks | No |

**Computation Formulas and Model**

**Level 1 — required minimum:**

Variables:

| Variable | Meaning | Units | Domain |
|---|---|---|---|
| `S` | Susceptible population | individuals | `S >= 0` |
| `E` | Exposed population (SEIR only) | individuals | `E >= 0` |
| `I` | Infectious population | individuals | `I >= 0` |
| `R` | Recovered population | individuals | `R >= 0` |
| `N` | Total population | individuals | `N > 0` (fixed, `N = 10000`) |
| `β` (beta) | Transmission rate | 1/day | `β > 0`, derived: `β = R₀ · γ` |
| `γ` (gamma) | Recovery rate | 1/day | `γ ∈ (0, 0.5]` |
| `σ` (sigma) | Incubation rate (SEIR only) | 1/day | `σ ∈ (0, 1.0]` |
| `R₀` | Basic reproduction number | dimensionless | `R₀ ∈ [0.5, 8.0]` |
| `v` | Vaccination coverage | fraction | `v ∈ [0, 0.9]` |
| `R_eff(t)` | Effective reproduction number | dimensionless | derived: `R_eff(t) = R₀ · S(t)/N` |

Relationships — **SIR model**:

```
dS/dt = -β·S·I/N
dI/dt = β·S·I/N - γ·I
dR/dt = γ·I
```

Relationships — **SEIR model**:

```
dS/dt = -β·S·I/N
dE/dt = β·S·I/N - σ·E
dI/dt = σ·E - γ·I
dR/dt = γ·I
```

Initial conditions:

```
S₀ = N·(1 - v) - I₀
I₀ = 1
E₀ = 0 (SEIR only)
R₀_init = N·v   (initial recovered = vaccinated fraction)
```

Where `v = vaccination / 100` (vaccination is provided as a percentage).

Derived outputs:

```
R_eff(t) = R₀ · S(t)/N
Peak infection: argmax_t I(t), max I(t)
Herd immunity threshold: HIT = 1 - 1/R₀
Final recovered: R(t_end)/N × 100%
```

Output properties (invariants):

| Property | Description |
|---|---|
| `S(t) + E(t) + I(t) + R(t) = N` for all t | Conservation of total population (SIR: S+I+R=N) |
| `S(t) >= 0, E(t) >= 0, I(t) >= 0, R(t) >= 0` | Non-negativity of all compartments |
| `S(t)` is monotonically non-increasing | Susceptible only decreases (no births/immigration) |
| `R(t)` is monotonically non-decreasing | Recovered only increases |
| `S(0) = N·(1-v) - 1, I(0) = 1, R(0) = N·v` | Initial conditions preserved |
| When `R₀ < 1`, `I(t)` decreases monotonically | No epidemic when R₀ < 1 |
| Peak of I(t) occurs when `R_eff = 1` | Mathematical property of the model |

Reference examples:

| # | Inputs | Expected output | Computational path | Source |
|---|---|---|---|---|
| 1 | SIR, R₀=3, γ=0.1, N=10000, v=0, S=9999, I=1, R=0 | dS/dt = -0.3·9999·1/10000 ≈ -0.300, dI/dt = 0.3·9999·1/10000 - 0.1·1 ≈ 0.200, dR/dt = 0.1 | SIR | Manual calculation |
| 2 | SEIR, R₀=3, γ=0.1, σ=0.2, N=10000, S=9999, E=0, I=1, R=0 | dS/dt ≈ -0.300, dE/dt ≈ 0.300, dI/dt = 0.2·0 - 0.1·1 = -0.1, dR/dt = 0.1 | SEIR | Manual calculation |
| 3 | SIR, R₀=3, γ=0.1, v=0.67, S=3299, I=1, R=6700 | dS/dt ≈ -0.0990, dI/dt ≈ -0.0011, dR/dt = 0.1. R_eff = 3·3299/10000 ≈ 0.99 < 1 → no epidemic | SIR, herd immunity | Manual calculation |
| 4 | SIR, R₀=1.5, γ=0.1, v=0 | HIT = 1 - 1/1.5 = 33.3% | HIT calculation | Formula |

**Level 2 — full formalization:**

- Complete mathematical formulation: Classical Kermack-McKendrick SIR model and its SEIR extension. Closed population (no births/deaths beyond disease), homogeneous mixing. IVP on [0, 300] days.
- Analytical properties: Disease-free equilibrium (S=N, I=0) is stable when R₀ < 1. Endemic threshold at R₀ = 1. Final size relation: `ln(S₀/S_∞) = R₀·(1 - S_∞/N)` (transcendental equation for final susceptibles). SEIR adds a delay but does not change the final size or R₀ threshold.
- Numerical method: MRT (Modified Rosenbrock Triple) from diff-grok, adaptive step, initial h=0.1. Suitable for both stiff and non-stiff ODEs.

Level 2 document location: in this specification.

> Level 2 need not be complete before implementation begins, but must be complete before the computational part is considered verified.

**Task input parameters:**

| Parameter | Type | Units | Domain | Description |
|---|---|---|---|---|
| `modelType` | `'SIR' \| 'SEIR'` | — | — | Model selection |
| `r0` | `number` | dimensionless | `[0.5, 8.0]` | Basic reproduction number |
| `gamma` | `number` | 1/day | `(0, 0.5]` | Recovery rate |
| `sigma` | `number` | 1/day | `(0, 1.0]` | Incubation rate (SEIR only) |
| `vaccination` | `number` | percent | `[0, 90]` | Initial vaccination coverage |

**Task output data:**

| Parameter | Type | Description |
|---|---|---|
| `t` | `Float64Array` | Time values (days) |
| `S` | `Float64Array` | Susceptible population over time |
| `E` | `Float64Array \| null` | Exposed population (SEIR) or null (SIR) |
| `I` | `Float64Array` | Infectious population over time |
| `R` | `Float64Array` | Recovered population over time |
| `rEff` | `Float64Array` | Effective reproduction number R_eff(t) |
| `peakDay` | `number` | Day of peak infection |
| `peakCount` | `number` | Peak number of infectious individuals |
| `finalRecoveredPct` | `number` | Final % of population recovered |
| `herdImmunityThreshold` | `number` | HIT = 1 - 1/R₀ |

**Computation implementation:**

| Step | Implementation method | Details | Documentation |
|---|---|---|---|
| 1 | External library | `diff-grok`, function `mrt(task: ODEs)`. MRT adaptive method. | [diff-grok](https://github.com/nicedoc/diff-grok) |
| 2 | Custom method | Compute R_eff(t) = R₀·S(t)/N for each time point | — |
| 3 | Custom method | Find peak: argmax I(t), max I(t) | — |
| 4 | Custom method | Summary: finalRecoveredPct = R(t_end)/N×100, HIT = 1 - 1/R₀ | — |

**Library call:**

```typescript
import { ODEs, mrt } from 'diff-grok';

// SIR: 3 equations; SEIR: 4 equations
const beta = r0 * gamma;

// SIR
const taskSIR: ODEs = {
  name: 'SIR',
  arg: { name: 't', start: 0, finish: 300, step: 0.1 },
  initial: [S0, I0, R0_init],
  func: (t, y, out) => {
    out[0] = -beta * y[0] * y[1] / N;           // dS/dt
    out[1] = beta * y[0] * y[1] / N - gamma * y[1]; // dI/dt
    out[2] = gamma * y[1];                        // dR/dt
  },
  tolerance: 1e-6,
  solutionColNames: ['S', 'I', 'R'],
};

// SEIR
const taskSEIR: ODEs = {
  name: 'SEIR',
  arg: { name: 't', start: 0, finish: 300, step: 0.1 },
  initial: [S0, E0, I0, R0_init],
  func: (t, y, out) => {
    out[0] = -beta * y[0] * y[2] / N;                // dS/dt
    out[1] = beta * y[0] * y[2] / N - sigma * y[1];  // dE/dt
    out[2] = sigma * y[1] - gamma * y[2];             // dI/dt
    out[3] = gamma * y[2];                             // dR/dt
  },
  tolerance: 1e-6,
  solutionColNames: ['S', 'E', 'I', 'R'],
};

const solution = mrt(modelType === 'SIR' ? taskSIR : taskSEIR);
```

**ODE right-hand side reference examples:**

| # | Model | R₀ | γ | σ | S | E | I | R | Expected dS/dt | Expected dE/dt | Expected dI/dt | Expected dR/dt | Derivation |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | SIR | 3.0 | 0.1 | — | 9999 | — | 1 | 0 | -0.02999 | — | 0.19999 | 0.1 | β=0.3; -0.3·9999·1/10000, 0.3·9999/10000-0.1 |
| 2 | SEIR | 3.0 | 0.1 | 0.2 | 9999 | 0 | 1 | 0 | -0.02999 | 0.02999 | -0.1 | 0.1 | β=0.3; dE=0.3·9999/10000-0, dI=0.2·0-0.1 |
| 3 | SIR | 3.0 | 0.1 | — | 5000 | — | 2000 | 3000 | -30.0 | — | 10.0 | 200.0 | β=0.3; -0.3·5000·2000/10000, 0.3·5000·2000/10000-0.1·2000 |
| 4 | SEIR | 2.0 | 0.2 | 0.5 | 8000 | 500 | 1000 | 500 | -32.0 | -218.0 | 150.0 | 200.0 | β=0.4; -0.4·8000·1000/10000, 0.4·8000·1000/10000-0.5·500, 0.5·500-0.2·1000 |
| 5 | SIR | 1.5 | 0.1 | — | 3300 | — | 1 | 6699 | -0.00495 | — | -0.09505 | 0.1 | β=0.15; R_eff=1.5·3300/10000=0.495<1, I decreases |

**Execution environment constraint:** Main thread (synchronous, < 100 ms). The core does not use Datagrok API — it depends only on `diff-grok` and custom methods.

---

### 1.2. Ports

#### Task ports: `task_primary`

| Port | Type | Interface / format | Description |
|---|---|---|---|
| Input | `EpidemicParams` | `{ modelType, r0, gamma, sigma, vaccination }` | 5 validated parameters |
| Output | `EpidemicSolution` | `{ t, S, E, I, R, rEff, peakDay, peakCount, finalRecoveredPct, herdImmunityThreshold }` | Solution arrays + summary |

#### Application-level ports

| Port | Used | Description |
|---|---|---|
| Progress | No | Computation is synchronous and fast (< 100 ms) |
| Cancellation | No | Not needed for fast synchronous computation |
| Data | No | No external data loading |

### 1.3. Adapters

| Adapter | Used | Implementation |
|---|---|---|
| UI adapter | Yes | Datagrok inputs (`ui.input.float`, `ui.input.toggle`), preset buttons, custom summary panel |
| Display adapter | Yes | Datagrok viewers (2 line charts, 1 scatter plot), custom summary panel (`comp_summary`) |
| Worker adapter | No | Computation is fast, main thread only |
| Progress adapter | No | Not needed |
| Data adapter | No | No external data loading |

### 1.4. Coordinator

The coordinator is the `sirSeirEpidemicApp()` function:

- **Input listening:** each input's `onValueChanged` calls `debouncedRun()`.
- **Reactivity management:** model type toggle controls visibility of σ slider and E(t) curve. R₀ slider recalculates β internally. Preset buttons batch-update multiple controls.
- **Validation trigger:** `runPrimary()` calls `validate(inputs)` before computation.
- **Control state during computations:** No blocking needed (computation < 100 ms).
- **Computation blocking:** `computationsBlocked` flag during preset application and reset to prevent multiple recomputations. Set `true` → update controls → set `false` → single `runPrimary()`.
- **Results to display:** Creates new `DG.DataFrame` with time series, assigns to view and viewers; updates summary panel labels.
- **Resource lifecycle:** `subs[]` collects subscriptions; cleaned up in `onViewRemoved`.

### 1.5. Independence Principle

Confirmed. Input controls (sliders, toggle, presets) and their reactivity (σ visibility, batch updates) function independently of the computational core. The core receives a ready `EpidemicParams` object with no knowledge of UI elements. The `validate()` function is pure, in `core.ts`.

---

## 2. Main View

| Field | Value |
|---|---|
| View type | `DG.TableView` |
| Description | Table view with docked line charts (epidemic dynamics, R_eff), scatter plot (phase portrait), left-panel form with controls, presets, and summary panel |

---

## 3. Controls (Inputs)

### 3.1. Primary Pipeline Controls

| ID | Label | Control type | Data type | Default | Min | Max | Step | Format | Nullable | Tooltip text | Group |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `ctrl_model_type` | Model | `ui.input.toggle` | `boolean` | `false` (SIR) | — | — | — | — | No | "Toggle between SIR (3 compartments) and SEIR (4 compartments, adds Exposed class with incubation period)" | Model |
| `ctrl_r0` | R₀ — basic reproduction number | `ui.input.float` | `number` | `3.0` | `0.5` | `8.0` | `0.1` | `0.0` | No | "Basic reproduction number — average number of people one infected person infects in a fully susceptible population. Influenza ≈ 1.5, COVID-19 ≈ 2.5–3.5, Measles ≈ 12–18." | Parameters |
| `ctrl_gamma` | γ — recovery rate (1/γ = infectious period in days) | `ui.input.float` | `number` | `0.1` | `0.01` | `0.5` | — | `0.00` | No | "1/γ is the average duration of the infectious period in days. For example, γ = 0.1 means a person is infectious for ~10 days on average." | Parameters |
| `ctrl_sigma` | σ — incubation rate (1/σ = incubation period in days) | `ui.input.float` | `number` | `0.2` | `0.05` | `1.0` | — | `0.00` | No | "1/σ is the average incubation period in days. During this time the person is infected (E) but not yet infectious to others." | Parameters |
| `ctrl_vaccination` | Initial vaccination coverage (%) | `ui.input.float` | `number` | `0` | `0` | `90` | — | `0` | No | "Initial fraction of immune individuals. Herd immunity threshold = 1 − 1/R₀. For R₀ = 3, at least 67% must be vaccinated to prevent an outbreak." | Parameters |

Note: `ctrl_sigma` is only visible/enabled when `ctrl_model_type` is SEIR (true).

### 3.2. Secondary Task Triggers

N/A — no secondary tasks.

### 3.3. Secondary Task Controls

N/A — no secondary tasks.

### 3.4. Other Buttons and Actions

| ID | Label / icon | Action | Tooltip text | Availability condition |
|---|---|---|---|---|
| `btn_influenza` | Influenza (R₀ ≈ 1.5) | Set R₀=1.5, γ=0.143 (7-day infectious period), σ=0.5 (2-day incubation) | "Influenza: R₀ ≈ 1.5, infectious period ~7 days, incubation ~2 days" | Always |
| `btn_covid` | COVID-19 (R₀ ≈ 3.0) | Set R₀=3.0, γ=0.1 (10-day infectious period), σ=0.2 (5-day incubation), switch to SEIR | "COVID-19: R₀ ≈ 3.0, infectious period ~10 days, incubation ~5 days. Switches to SEIR model." | Always |
| `btn_measles` | Measles (R₀ ≈ 15) | Set R₀=8.0 (clamped to max), γ=0.125 (8-day infectious period), σ=0.1 (10-day incubation), switch to SEIR | "Measles: R₀ ≈ 15 (clamped to slider max 8.0), infectious period ~8 days, incubation ~10 days. Switches to SEIR model." | Always |
| `btn_reset` | `ui.iconFA('undo')` | Reset all controls to defaults | "Reset all parameters to default values" | Always |

### 3.5. Custom UI Components

| Component ID | Brief description | Role | UI component specification |
|---|---|---|---|
| `comp_summary` | Summary panel showing epidemic key metrics. Built with `ui.divV`/`ui.label`. Displays: peak infection (day + count), final % recovered, herd immunity threshold, current R₀. Values updated via `textContent`. Styled via `.sir-seir-app-summary-panel`. | Display | Inline |

**`comp_summary` structure:**

```
ui.h2('Epidemic Summary')
  "Peak infection: day {peakDay}, {peakCount} cases"
  "Final recovered: {finalRecoveredPct}% of population"
  "Herd immunity threshold: {HIT}% of population"
  "Basic reproduction number R₀ = {r0}"
```

---

## 4. Result Display Elements

### 4.1. Primary Pipeline Display Elements

| ID | Type | Associated output data | Docking location |
|---|---|---|---|
| `view_epidemic` | Datagrok viewer `line chart` | `Time`, `S`, `E` (SEIR), `I`, `R` | Right area, top |
| `view_reff` | Datagrok viewer `line chart` | `Time`, `R_eff` | Right area, middle |
| `view_phase` | Datagrok viewer `scatter plot` | `S`, `I` | Right area, bottom |
| `comp_summary` | Custom HTMLElement | Summary statistics | Left panel, below controls |
| `view_table` | `DG.TableView` grid | All columns | Center (default) |

**Details for `view_epidemic`:**

| Property | Value |
|---|---|
| Title | Epidemic Dynamics — Population Over Time |
| X axis | `Time (days since first case)` |
| Y axis | `Number of individuals` |
| Series | S — blue line ("S — Susceptible (not yet infected)"), E — orange line ("E — Exposed (incubating, SEIR only)", hidden when SIR), I — red line ("I — Infectious (can spread)"), R — green line ("R — Recovered (immune)") |
| Special | Area fill under I(t) in light red to emphasize infection peak. Vertical dashed line at peak of I(t) with label "Peak: day N, I = X". |

**Details for `view_reff`:**

| Property | Value |
|---|---|
| Title | Effective Reproduction Number R_eff(t) |
| X axis | `Time (days since first case)` |
| Y axis | `R_eff` |
| Series | R_eff — single line |
| Special | Horizontal reference line at R_eff = 1 labeled "Epidemic threshold (R_eff = 1)". |

**Details for `view_phase`:**

| Property | Value |
|---|---|
| Title | Phase Portrait — Susceptible vs Infectious |
| X axis | `S — Susceptible` |
| Y axis | `I — Infectious` |
| Special | Trajectory from initial point through peak back to I=0. Peak point marked with distinct marker. |

### 4.2. Secondary Task Display Elements

N/A — no secondary tasks.

---

## 5. Layout and UI Element Placement

### 5.1. Control Placement

| Area | Content (control IDs) |
|---|---|
| Left panel | `ui.form` with grouped controls: Model group (`ctrl_model_type`), Parameters group (`ctrl_r0`, `ctrl_gamma`, `ctrl_sigma`, `ctrl_vaccination`), Disease Presets group (`btn_influenza`, `btn_covid`, `btn_measles`), Epidemic Summary (`comp_summary`) |
| Ribbon (Actions group) | `btn_reset` |
| Main area | `view_table` (grid) |

**Structure of left panel form:**

```
ui.h2('Model')
  ctrl_model_type (SIR / SEIR toggle)

ui.h2('Parameters')
  ctrl_r0
  ctrl_gamma
  ctrl_sigma  (visible only when SEIR)
  ctrl_vaccination

ui.h2('Disease Presets')
  btn_influenza
  btn_covid
  btn_measles

ui.h2('Epidemic Summary')
  comp_summary
```

### 5.2. Display Element Placement

| Element ID | Docking area | Position / ratio |
|---|---|---|
| `form` (left panel) | `DG.DOCK_TYPE.LEFT` relative to root | ratio `0.6` |
| `view_epidemic` | `DG.DOCK_TYPE.RIGHT` relative to root | ratio `0.7` |
| `view_reff` | `DG.DOCK_TYPE.DOWN` relative to `view_epidemic` | ratio `0.35` |
| `view_phase` | `DG.DOCK_TYPE.DOWN` relative to `view_reff` | ratio `0.5` |
| `view_table` (grid) | Default position (center) | — |

### 5.3. Styles

CSS file: `css/sir-seir.css`. Import: `import '../../css/sir-seir.css'`. All custom classes use the `sir-seir-app-` prefix for isolation.

**Static styles:**

| Element | CSS class(es) | Description |
|---|---|---|
| Summary panel | `.sir-seir-app-summary-panel` | Background, padding, border-radius, gap |
| Summary section header | `.sir-seir-app-summary-panel > label` | Bold, grey color |
| Summary value | `.sir-seir-app-summary-value` | Monospace font, right-aligned |
| Preset button group | `.sir-seir-app-preset-group` | Flex row, gap between buttons |

**Dynamic styles:**

| Element | CSS class(es) | Condition | Description |
|---|---|---|---|
| σ slider container | `.sir-seir-app-sigma-hidden` | When model is SIR | `display: none` to hide σ control |

---

## 6. User Feedback

### 6.1. Control Tooltips

| Control ID | Tooltip text | Mechanism |
|---|---|---|
| `ctrl_model_type` | Toggle between SIR (3 compartments) and SEIR (4 compartments, adds Exposed class with incubation period) | `tooltipText` property |
| `ctrl_r0` | Basic reproduction number — average number of people one infected person infects in a fully susceptible population. Influenza ≈ 1.5, COVID-19 ≈ 2.5–3.5, Measles ≈ 12–18. | `tooltipText` property |
| `ctrl_gamma` | 1/γ is the average duration of the infectious period in days. For example, γ = 0.1 means a person is infectious for ~10 days on average. | `tooltipText` property |
| `ctrl_sigma` | 1/σ is the average incubation period in days. During this time the person is infected (E) but not yet infectious to others. | `tooltipText` property |
| `ctrl_vaccination` | Initial fraction of immune individuals. Herd immunity threshold = 1 − 1/R₀. For R₀ = 3, at least 67% must be vaccinated to prevent an outbreak. | `tooltipText` property |
| `btn_influenza` | Influenza: R₀ ≈ 1.5, infectious period ~7 days, incubation ~2 days | Button tooltip |
| `btn_covid` | COVID-19: R₀ ≈ 3.0, infectious period ~10 days, incubation ~5 days. Switches to SEIR model. | Button tooltip |
| `btn_measles` | Measles: R₀ ≈ 15 (clamped to slider max 8.0), infectious period ~8 days, incubation ~10 days. Switches to SEIR model. | Button tooltip |
| `btn_reset` | Reset all parameters to default values | `ui.iconFA` third argument |

### 6.2. Validators as Feedback

| Input ID | Validation source | Description |
|---|---|---|
| `ctrl_r0` … `ctrl_vaccination` | `validate()` from core | Each input has `addValidator()` calling `validate(getInputs())` and returning the error for its own ID, or `null`. Displays inline hint on invalid value. |

### 6.3. Progress Bar

| Task | Progress bar | Type | Cancellation support |
|---|---|---|---|
| `task_primary` | No | — | — |

---

## 7. Validation

### 7.1. Primary Pipeline Validation

#### Complex Validation Rules

| Rule ID | Condition (invalid) | Affected inputs (ID) | Error message |
|---|---|---|---|
| `val_01` | `r0 < 0.5 \|\| r0 > 8.0` | `ctrl_r0` | "R₀ must be between 0.5 and 8.0" |
| `val_02` | `gamma <= 0 \|\| gamma > 0.5` | `ctrl_gamma` | "Recovery rate must be in (0, 0.5]" |
| `val_03` | `sigma <= 0 \|\| sigma > 1.0` (SEIR only) | `ctrl_sigma` | "Incubation rate must be in (0, 1.0]" |
| `val_04` | `vaccination < 0 \|\| vaccination > 90` | `ctrl_vaccination` | "Vaccination coverage must be between 0% and 90%" |
| `val_05` | `vaccination / 100 >= 1 - 1/N` (all vaccinated, no susceptibles left for I₀=1) | `ctrl_vaccination` | "Vaccination coverage too high — no susceptible individuals remain" |

#### Validation Order

```
1. val_01 (R₀ range)
2. val_02 (γ range)
3. val_03 (σ range — skipped when SIR)
4. val_04 (vaccination range)
5. val_05 (vaccination feasibility — checked only if val_04 passed)
```

#### Returned Map Format

```
Map<InputId, string>
```

Where `InputId = 'ctrl_r0' | 'ctrl_gamma' | 'ctrl_sigma' | 'ctrl_vaccination'`.

### 7.2. Secondary Task Validation

N/A — no secondary tasks.

---

## 8. Main Pipeline

### 8.1. Primary Pipeline

| Step | Description |
|---|---|
| Parameter input | User sets values through controls → `getInputs()` converts to `EpidemicParams` |
| Validation | `validate(inputs)` returns `Map<InputId, string>`. On errors: mark invalid inputs, call `clearResults()` |
| Computation | `solve(inputs)` — synchronous ODE solution via `mrt()` + derived quantities |
| Result display | `updateDataFrame(result)` creates new DataFrame, assigns to view and viewers; `updateSummaryPanel(result)` updates labels |

Reactive trigger: `onValueChanged` with debounce 50 ms.

Error behavior on validation failure: clear results (empty DataFrame).

Error behavior on computation failure: clear results + `grok.shell.error(msg)`.

### 8.2. Secondary Pipelines

N/A — no secondary tasks.

### 8.3. Common Pipeline Aspects

#### Control Behavior During Computations

| Pipeline / task | Controls blocked | Which controls |
|---|---|---|
| `task_primary` | No | — (computation < 100 ms) |

#### Computation Error Handling

| Pipeline / task | Strategy | Notification method |
|---|---|---|
| `task_primary` | Clear results | `grok.shell.error` |

### 8.4. Computation Blocking and Batch Input Updates

| Scenario | Source | Target controls (ID) | Blocked pipelines |
|---|---|---|---|
| Disease preset application | `btn_influenza` / `btn_covid` / `btn_measles` | `ctrl_r0`, `ctrl_gamma`, `ctrl_sigma`, `ctrl_model_type` | Primary (blocked) |
| Reset to defaults | `btn_reset` | All `ctrl_*` | Primary (blocked) |
| Format setting on init | Coordinator | All `ctrl_*` | Primary (blocked) |

Reactivity mode during batch update:

| Scenario | Reactivity mode |
|---|---|
| All above | `computationsBlocked = true` → writes → `computationsBlocked = false` → single `runPrimary()` |

---

## 9. Reactivity and Dependencies Between Inputs

### 9.1. Dependency Graph

| Source (input ID) | Target (input IDs) | Reaction type | Logic |
|---|---|---|---|
| `ctrl_model_type` | `ctrl_sigma` | Availability | σ slider visible/enabled only when SEIR is selected |
| `ctrl_model_type` | `view_epidemic` | Display update | Show/hide E(t) series |
| Any `ctrl_*` | All viewers + `comp_summary` | Reactive update | Rerun `task_primary`, update all outputs |

### 9.2. Debounce / Throttle

| Input ID | Strategy | Interval (ms) |
|---|---|---|
| All `ctrl_*` | debounce | 50 |

---

## 10. Data Lifecycle

### 10.1. Data Input

Primary method: manual input via `ui.form` controls.

Initial state: all controls initialized with defaults (SIR, R₀=3, γ=0.1, σ=0.2, v=0%). `task_primary` runs automatically on initialization.

### 10.2. Loading from Resources

N/A — no external data loading.

### 10.3. Results Table Lifecycle

Update strategy: **DataFrame replacement** — on each recomputation a new `DG.DataFrame` is created and assigned to `view.dataFrame` and all viewers.

```
1. Application initialization
   → solve with defaults → new DG.DataFrame with columns [Time, S, I, R, R_eff] (SIR)
     or [Time, S, E, I, R, R_eff] (SEIR)
   → DataFrame added to TableView

2. task_primary completion
   → new DataFrame created from solution arrays
   → assigned to view.dataFrame, all viewer.dataFrame
   → all viewers redrawn reactively

3. Model type toggle (SIR ↔ SEIR)
   → DataFrame columns change (E column added/removed)
   → new DataFrame created with appropriate columns

4. Preset button press
   → batch update controls → single task_primary → new DataFrame

5. btn_reset press
   → batch update controls → single task_primary → new DataFrame
```

---

## 11. Error Handling Beyond Computations

| Error type | Strategy | Notification method |
|---|---|---|
| ODE solver error | Clear results, show message | `grok.shell.error` |
| ODE solver returns NaN/Inf | Clear results, show message | `grok.shell.error("Numerical instability detected")` |

Data loading, network, and worker errors are N/A (no external data, no workers).

---

## 12. Subscriptions and Resource Management

### 12.1. Event Subscriptions

| Subscription | Event | Cleanup mechanism |
|---|---|---|
| View removal | `grok.events.onViewRemoved` | `sub.unsubscribe()` in cleanup handler |
| Input value changes | `onValueChanged` on each `ctrl_*` | `sub.unsubscribe()` via `subs[]` |

All subscriptions collected in `subs[]` and unsubscribed when the view is removed.

### 12.2. Worker Termination

N/A — no web workers.

---

## 13. Application Closure

On view close (`grok.events.onViewRemoved`), the coordinator performs:

- [x] All event subscriptions unsubscribed (`subs[]` iteration)
- [x] Pending debounce timer cleared (`clearTimeout(debounceTimer)`)
- [x] No web workers to terminate
- [x] No open secondary dialogs

Closure handler: `grok.events.onViewRemoved.subscribe(...)`

---

## 14. Accessibility and UX

### 14.1. Keyboard Shortcuts

N/A — no custom keyboard shortcuts.

### 14.2. Context Menus

N/A — no custom context menus.

### 14.3. Undo / Redo

Not supported. `btn_reset` resets all controls to default values as the only rollback mechanism.

---

## 15. Testing

### 15.1. Computational Part (Core)

| File | Categories | Test count | Description |
|---|---|---|---|
| `src/tests/sir-seir-math-tests.ts` | Math: SIR func, Math: SEIR func, Math: SIR solve properties, Math: SEIR solve properties, Math: R_eff, Math: Summary stats | ~20 | ODE RHS verification, solve properties, R_eff calculation, summary statistics |
| `src/tests/sir-seir-validation-tests.ts` | API: Validation | ~14 | Validation rules for all inputs, defaults, boundaries, edge cases |

Tests are run via `grok test` (entry point: `src/package-test.ts`).

### 15.2. Inputs

| Category | Coverage | Description |
|---|---|---|
| Boundary values | `v_bnd_r0_low` (R₀=0.5), `v_bnd_gamma_max` (γ=0.5), `v_bnd_sigma_max` (σ=1.0) | Allowed boundary values |
| Invalid values | `v_01`–`v_05` (10 tests) | Out-of-range and edge cases for each parameter |
| Model-dependent | `v_sigma_sir` | σ validation skipped when SIR mode |
| Multiple simultaneous errors | `v_multi` | r0=0, gamma=0, vaccination=100 → ≥ 3 errors |
| Valid defaults | `v_def` | All defaults pass validation |

### 15.3. Mathematical Verification

#### Level 1 Verification (required)

**Formula/equation verification:**

| Test category | Test count | What is verified | Reference source |
|---|---|---|---|
| Math: SIR func | 3 | ODE RHS: concrete (β, γ, S, I, R) → expected (dS/dt, dI/dt, dR/dt) | Manual calculation |
| Math: SEIR func | 3 | ODE RHS: concrete (β, γ, σ, S, E, I, R) → expected (dS/dt, dE/dt, dI/dt, dR/dt) | Manual calculation |
| Math: R_eff | 2 | R_eff = R₀·S/N for specific values | Manual calculation |
| Math: Summary stats | 2 | HIT = 1-1/R₀, final recovered % | Manual calculation |

**Output property verification:**

| Test category | Test count | Properties verified |
|---|---|---|
| Math: SIR solve properties | 4 | Non-negativity, initial conditions, population conservation (S+I+R=N), S monotonically non-increasing |
| Math: SEIR solve properties | 4 | Non-negativity, initial conditions, population conservation (S+E+I+R=N), S monotonically non-increasing |
| Math: R₀ < 1 behavior | 1 | I(t) decreases monotonically when R₀ < 1 |
| Math: Herd immunity | 1 | Near-zero epidemic when vaccination ≥ HIT |

#### Level 2 Verification (for full formalization)

**Numerical method verification:**

| Test category | Test count | Reference problems | Tolerance | Source |
|---|---|---|---|---|
| Math: MRT solver | 2 | Non-stiff 1D (`dy/dt = 4·exp(0.8t) − 0.5y`), Stiff 1D (`dy/dt = −1000y + 3000 − 2000·exp(−t)`) | Max error < 0.1 | Chapra & Canale textbook |

**Convergence verification:**

| Status | Description |
|---|---|
| Not yet covered | Could verify that reducing tolerance decreases solution discrepancy |

**Asymptotic/equilibrium behavior:**

| Status | Description |
|---|---|
| Implemented | For R₀ > 1: verify I(t_end) ≈ 0 (epidemic ends). For R₀ < 1: verify I(t) → 0 monotonically. |

#### Validation Tests

| Test ID | Rule | Input data | Expected result |
|---|---|---|---|
| `v_01a` | val_01 | `r0 = 0.4` | Error on `ctrl_r0` |
| `v_01b` | val_01 | `r0 = 8.5` | Error on `ctrl_r0` |
| `v_02a` | val_02 | `gamma = 0` | Error on `ctrl_gamma` |
| `v_02b` | val_02 | `gamma = -0.1` | Error on `ctrl_gamma` |
| `v_02c` | val_02 | `gamma = 0.6` | Error on `ctrl_gamma` |
| `v_03a` | val_03 | `sigma = 0` (SEIR) | Error on `ctrl_sigma` |
| `v_03b` | val_03 | `sigma = 1.5` (SEIR) | Error on `ctrl_sigma` |
| `v_03c` | val_03 | `sigma = 0` (SIR) | No error (skipped) |
| `v_04a` | val_04 | `vaccination = -5` | Error on `ctrl_vaccination` |
| `v_04b` | val_04 | `vaccination = 95` | Error on `ctrl_vaccination` |
| `v_05a` | val_05 | `vaccination = 90`, N=10000 | Error on `ctrl_vaccination` (no susceptibles for I₀=1) |
| `v_def` | all defaults | Defaults | `errors.size = 0` |
| `v_bnd_r0` | valid boundary | `r0 = 0.5` | No error |
| `v_multi` | multiple | `r0=0, gamma=0, vaccination=100` | ≥ 3 errors |

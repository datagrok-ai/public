# Application Specification: Circular Restricted Three-Body Problem (CR3BP)

<!-- TOC for partial reading -->
<!-- Section 1: General Architecture — lines 15–220 -->
<!-- Section 2: Main View — lines 221–229 -->
<!-- Section 3: Controls — lines 230–310 -->
<!-- Section 4: Display Elements — lines 311–380 -->
<!-- Section 5: Layout — lines 381–445 -->
<!-- Section 6: User Feedback — lines 446–490 -->
<!-- Section 7: Validation — lines 491–530 -->
<!-- Section 8: Pipeline — lines 531–610 -->
<!-- Section 9: Reactivity — lines 611–635 -->
<!-- Sections 10–15: Data, Errors, Resources, Closure, UX, Testing — lines 636–850 -->

## 1. General Architecture

### 1.0. General Information

| Field | Value |
|---|---|
| Application name | Circular Restricted Three-Body Problem |
| Package | ThreeBodyProblem |
| Entry function | `threeBodyProblemApp()` |
| Brief description | Interactive simulation of the Circular Restricted Three-Body Problem (CR3BP) in the co-rotating frame. Visualizes orbits, Lagrange points, Jacobi constant conservation, and zero-velocity (Hill) surfaces for user-defined mass ratios and initial conditions. |
| Main view type | `DG.TableView` |

### 1.1. Core

#### Task List

| Task ID | Name | Pipeline type | Trigger | Synchronicity | Execution environment | Parallelization |
|---|---|---|---|---|---|---|
| `task_primary` | CR3BP Integration & Analysis | Primary (reactive) | Any input change | Sync | Main thread | No |

#### Dependencies Between Tasks

```
task_primary — independent (single task)
```

#### Task: `task_primary`

**General characteristics:**

| Field | Value |
|---|---|
| Name | CR3BP Integration & Analysis |
| Description | Numerically integrates the CR3BP equations of motion, computes all 5 Lagrange points, the Jacobi constant along the trajectory, and a 2D grid of the effective potential for zero-velocity curve visualization |
| Dependency on other tasks | No |

**Computation Formulas and Model**

**Level 1 — required minimum:**

Variables:

| Variable | Meaning | Units | Domain |
|---|---|---|---|
| `x` | Position along primary axis (synodic frame) | dimensionless | `x ∈ ℝ` |
| `y` | Position perpendicular to axis (synodic frame) | dimensionless | `y ∈ ℝ` |
| `vx` | Velocity in x direction (synodic frame) | dimensionless | `vx ∈ ℝ` |
| `vy` | Velocity in y direction (synodic frame) | dimensionless | `vy ∈ ℝ` |
| `μ` | Mass parameter m₂/(m₁+m₂) | dimensionless | `μ ∈ (0, 0.5]` |
| `x₀` | Initial x position | dimensionless | `x₀ ∈ ℝ` |
| `y₀` | Initial y position | dimensionless | `y₀ ∈ ℝ` |
| `vx₀` | Initial x velocity | dimensionless | `vx₀ ∈ ℝ` |
| `vy₀` | Initial y velocity | dimensionless | `vy₀ ∈ ℝ` |
| `T` | Integration time | dimensionless time units | `T > 0` |
| `r₁` | Distance to body 1 (mass 1−μ, at (−μ, 0)) | dimensionless | `r₁ > 0` |
| `r₂` | Distance to body 2 (mass μ, at (1−μ, 0)) | dimensionless | `r₂ > 0` |

Relationships (equations of motion in the synodic frame):

```
r₁ = √((x + μ)² + y²)
r₂ = √((x − 1 + μ)² + y²)

dx/dt  = vx
dy/dt  = vy
dvx/dt = 2·vy + x − (1−μ)·(x+μ)/r₁³ − μ·(x−1+μ)/r₂³
dvy/dt = −2·vx + y − (1−μ)·y/r₁³ − μ·y/r₂³
```

Jacobi constant (integral of motion):

```
U(x,y) = 0.5·(x² + y²) + (1−μ)/r₁ + μ/r₂
C_J = 2·U(x,y) − (vx² + vy²)
```

Lagrange points (equilibrium points where all accelerations = 0):

- **L1, L2, L3** (collinear, on x-axis, y=0): Found by solving the quintic equation for x with y=0, vx=0, vy=0. Numerically solved via Newton's method on `f(x) = x − (1−μ)·(x+μ)/|x+μ|³ − μ·(x−1+μ)/|x−1+μ|³ = 0`.
  - L1: between the two bodies, `x ∈ (−μ, 1−μ)`
  - L2: beyond the smaller body, `x > 1−μ`
  - L3: beyond the larger body, `x < −μ`

- **L4, L5** (triangular): Form equilateral triangles with the two bodies.
  - L4: `(0.5 − μ, √3/2)`
  - L5: `(0.5 − μ, −√3/2)`

Zero-velocity curves: The boundary where `C_J = 2·U(x,y)` (i.e., velocity = 0). For a given `C_J` value, the forbidden region is where `2·U(x,y) < C_J` (kinetic energy would be negative).

Output properties (invariants):

| Property | Description |
|---|---|
| `C_J(t) = const` | Jacobi constant is conserved — drift indicates numerical error |
| `C_J(0)` determines forbidden regions | Zero-velocity curves partition space into accessible/forbidden regions |
| L4, L5 stable when `μ < 0.0385…` | Routh's criterion for triangular point stability |
| L1, L2, L3 always unstable | Collinear points are saddle points of the effective potential |

Reference examples:

| # | Inputs | Expected output | Computational path | Source |
|---|---|---|---|---|
| 1 | `μ=0.01215, x=0.8, y=0, vx=0, vy=0` | `r₁=√((0.8+0.01215)²)=0.81215`, `r₂=√((0.8−0.98785)²)=0.18785`, `dvx/dt=2·0+0.8−(0.98785·0.81215/0.81215³)−(0.01215·(−0.18785)/0.18785³)` | ODE RHS evaluation | Manual calculation |
| 2 | `μ=0.01215` | `L4=(0.48785, 0.86603)`, `L5=(0.48785, −0.86603)` | Triangular Lagrange points | Analytical: `(0.5−μ, ±√3/2)` |
| 3 | `μ=0.01215, x=0.8369, y=0, vx=0, vy=0` | `C_J ≈ 3.188` (near L1 value) | Jacobi constant | Literature (Earth-Moon L1) |

**Level 2 — full formalization:**

- Complete mathematical formulation: CR3BP is a Hamiltonian system in the rotating frame with 2 degrees of freedom. The synodic frame rotates with angular velocity 1. Distances normalized so that the separation between primaries = 1, total mass = 1.
- Analytical properties: 5 equilibrium points (3 collinear, 2 triangular). Collinear points are unstable (index-1 saddles). Triangular points are linearly stable for `μ < μ_Routh ≈ 0.03852`. The Jacobi constant is the only integral of motion.
- Numerical method: MRT (Modified Rosenbrock Triple) — an adaptive implicit method suitable for the mildly stiff CR3BP dynamics near close approaches. Initial step `h=0.001`, tolerance `1e-8`.

Level 2 document location: in this specification.

**Task input parameters:**

| Parameter | Type | Units | Domain | Description |
|---|---|---|---|---|
| `mu` | `number` | dimensionless | `(0, 0.5]` | Mass parameter |
| `x0` | `number` | dimensionless | `ℝ` | Initial x position |
| `y0` | `number` | dimensionless | `ℝ` | Initial y position |
| `vx0` | `number` | dimensionless | `ℝ` | Initial x velocity |
| `vy0` | `number` | dimensionless | `ℝ` | Initial y velocity |
| `T` | `number` | dimensionless time | `> 0` | Integration time |

**Task output data:**

| Parameter | Type | Description |
|---|---|---|
| `t` | `Float64Array` | Time values |
| `x` | `Float64Array` | x position trajectory |
| `y` | `Float64Array` | y position trajectory |
| `vx` | `Float64Array` | x velocity trajectory |
| `vy` | `Float64Array` | y velocity trajectory |
| `vmag` | `Float64Array` | Velocity magnitude at each step: `√(vx² + vy²)` |
| `cj` | `Float64Array` | Jacobi constant at each step |
| `lagrangePoints` | `{name: string, x: number, y: number}[]` | All 5 Lagrange points (L1–L5) |
| `zvcGrid` | `{xArr: Float64Array, yArr: Float64Array, U: Float64Array, cj0: number}` | Grid of effective potential values for ZVC plotting |

**Computation implementation:**

| Step | Implementation method | Details | Documentation |
|---|---|---|---|
| 1 | External library | `diff-grok`, function `mrt(task: ODEs)`. MRT adaptive method for the 4D CR3BP system. | [diff-grok](https://github.com/datagrok-ai/diff-grok) |
| 2 | Custom method | Lagrange point computation: L4/L5 analytical, L1/L2/L3 via Newton's method on the collinear equilibrium equation | — |
| 3 | Custom method | Jacobi constant: `C_J(t) = 2·U(x,y) − (vx² + vy²)` for each time step | — |
| 4 | Custom method | Zero-velocity curve grid: evaluate `2·U(x,y)` on a 200×200 grid over `[−1.5, 1.5] × [−1.5, 1.5]` | — |

**Library call:**

```typescript
import { ODEs, mrt } from 'diff-grok';

const task: ODEs = {
  name: 'CR3BP',
  arg: { name: 't', start: 0, finish: T, step: 0.001 },
  initial: [x0, y0, vx0, vy0],
  func: (t, y, out) => {
    const r1 = Math.sqrt((y[0] + mu) ** 2 + y[1] ** 2);
    const r2 = Math.sqrt((y[0] - 1 + mu) ** 2 + y[1] ** 2);
    const r1_3 = r1 ** 3;
    const r2_3 = r2 ** 3;
    out[0] = y[2]; // dx/dt = vx
    out[1] = y[3]; // dy/dt = vy
    out[2] = 2 * y[3] + y[0] - (1 - mu) * (y[0] + mu) / r1_3 - mu * (y[0] - 1 + mu) / r2_3;
    out[3] = -2 * y[2] + y[1] - (1 - mu) * y[1] / r1_3 - mu * y[1] / r2_3;
  },
  tolerance: 1e-8,
  solutionColNames: ['x', 'y', 'vx', 'vy'],
};
const solution = mrt(task);
```

**ODE right-hand side reference examples:**

| # | μ | x | y | vx | vy | Expected dx/dt | Expected dy/dt | Expected dvx/dt | Expected dvy/dt | Derivation |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.01215 | 0.8 | 0 | 0 | 0 | 0 | 0 | `2·0 + 0.8 − 0.98785·0.81215/0.81215³ − 0.01215·(−0.18785)/0.18785³` | `−2·0 + 0 − 0` | Substitution into EOM |
| 2 | 0.01215 | 0.48785 | 0.86603 | 0 | 0 | 0 | 0 | ≈ 0.48785 (L4 x-accel component) | ≈ 0.86603 (L4 y-accel component) | At L4, net acceleration = Coriolis only when v=0; true equilibrium requires vx=vy=0 |
| 3 | 0.5 | 0 | 0.86603 | 0 | 0 | 0 | 0 | ≈ 0 | ≈ 0.86603 | Equal mass case, L4 point |

**Execution environment constraint:** main thread (synchronous). The core does not use Datagrok API — depends only on `diff-grok` and custom methods.

---

### 1.2. Ports

#### Task ports: `task_primary`

| Port | Type | Interface / format | Description |
|---|---|---|---|
| Input | `CR3BPParams` | `{ mu, x0, y0, vx0, vy0, T }` | 6 validated parameters |
| Output | `CR3BPSolution` | `{ t, x, y, vx, vy, vmag, cj, lagrangePoints, zvcGrid }` | Trajectory + derived data |

#### Application-level ports

| Port | Used | Description |
|---|---|---|
| Progress | No | Computation is synchronous and fast (< 500ms) |
| Cancellation | No | No long-running async tasks |
| Data | No | No external data loading |

### 1.3. Adapters

| Adapter | Used | Implementation |
|---|---|---|
| UI adapter | Yes | Datagrok inputs (`ui.input.float` × 5, `ui.input.float` × 1), preset buttons, click-to-launch on scatter plot |
| Display adapter | Yes | Datagrok viewers (scatter plot × 2, line chart × 1), custom Lagrange point markers via `onAfterDrawScene` |
| Worker adapter | No | — |
| Progress adapter | No | — |
| Data adapter | No | — |

### 1.4. Coordinator

The coordinator is the `threeBodyProblemApp()` function:

- **Input listening:** each control calls `debouncedRun()` from `onValueChanged`. Scatter plot click handler sets `x0, y0` controls (and optionally `vx0, vy0` on drag).
- **Reactivity management:** all controls trigger the same primary pipeline; preset buttons perform batch updates.
- **Validation trigger:** `runPrimary()` calls `validate(inputs)` before computation.
- **Control state during computations:** no blocking needed (< 500ms synchronous).
- **Computation blocking:** `computationsBlocked = true` during preset button batch writes, then single `runPrimary()`.
- **Results to display:** `updateDataFrame()` creates new `DG.DataFrame` and assigns to view and all viewers; Lagrange points and ZVC data rendered as overlays or separate DataFrame.
- **Resource lifecycle:** `subs[]` collects subscriptions; cleaned up in `onViewRemoved`.

### 1.5. Independence Principle

Confirmed. Input controls and their reactivity (debounce, presets) function independently of the computational core. The core receives a `CR3BPParams` object and returns a `CR3BPSolution`. The `validate()` function is pure and in `core.ts`. Click-to-launch interaction only manipulates control values — it does not call the core directly.

---

## 2. Main View

| Field | Value |
|---|---|
| View type | `DG.TableView` |
| Description | Table view with docked scatter plots (orbit, ZVC), line chart (Jacobi constant), and left-panel form with controls, presets, and Lagrange point info |

---

## 3. Controls (Inputs)

### 3.1. Primary Pipeline Controls

| ID | Label | Control type | Data type | Default | Min | Max | Format | Nullable | Tooltip text | Group |
|---|---|---|---|---|---|---|---|---|---|---|
| `ctrl_mu` | μ — mass ratio m₂/(m₁+m₂) | `ui.input.float` | `number` | `0.01215` | `0.001` | `0.5` | `0.00000` | No | Ratio of the smaller body's mass to the total. For Earth–Moon, μ ≈ 0.012. Increasing μ enlarges the gravitational sphere of influence of the second body — the topology of forbidden regions changes. | System Parameters |
| `ctrl_x0` | Initial position x₀ | `ui.input.float` | `number` | `0.994` | `-1.5` | `1.5` | `0.000` | No | Starting x coordinate of the test body in the rotating frame. Can also be set by clicking on the orbit plot. | Initial Conditions |
| `ctrl_y0` | Initial position y₀ | `ui.input.float` | `number` | `0.0` | `-1.5` | `1.5` | `0.000` | No | Starting y coordinate of the test body in the rotating frame. Can also be set by clicking on the orbit plot. | Initial Conditions |
| `ctrl_vx0` | Initial velocity vx₀ (rotating frame) | `ui.input.float` | `number` | `0.0` | `-2.0` | `2.0` | `0.000` | No | Initial velocity in the x direction. Can be set by dragging from the click point on the orbit plot. | Initial Conditions |
| `ctrl_vy0` | Initial velocity vy₀ (rotating frame) | `ui.input.float` | `number` | `-2.0016` | `-2.0` | `2.0` | `0.000` | No | Initial velocity in the y direction. Can be set by dragging from the click point on the orbit plot. | Initial Conditions |
| `ctrl_T` | Integration time (orbital periods) | `ui.input.float` | `number` | `20.0` | `1.0` | `100.0` | `0.0` | No | Total integration time in dimensionless units. One lunar orbital period ≈ 2π ≈ 6.28 time units. | Integration |

### 3.2. Secondary Task Triggers

N/A — no secondary tasks.

### 3.3. Secondary Task Controls

N/A — no secondary tasks.

### 3.4. Other Buttons and Actions

| ID | Label / icon | Action | Tooltip text | Availability condition |
|---|---|---|---|---|
| `btn_preset_earth_moon` | Earth–Moon (μ = 0.01215) | Set `ctrl_mu = 0.01215` | Set mass parameter to the Earth–Moon system value | Always |
| `btn_preset_sun_jupiter` | Sun–Jupiter (μ = 0.000953) | Set `ctrl_mu = 0.000953` | Set mass parameter to the Sun–Jupiter system value | Always |
| `btn_preset_pluto_charon` | Pluto–Charon (μ = 0.1) | Set `ctrl_mu = 0.1` | Set mass parameter to the Pluto–Charon system value | Always |
| `btn_preset_halo` | Earth–Moon L1 Halo Orbit | Batch set: `ctrl_mu=0.01215, ctrl_x0=0.8369, ctrl_y0=0, ctrl_vx0=0, ctrl_vy0=-0.0559, ctrl_T=6.19` | A periodic orbit around L1 — used by real space missions (e.g., JWST). Requires precise initial conditions. | Always |
| `btn_preset_freereturn` | Free-Return Trajectory | Batch set: `ctrl_mu=0.01215, ctrl_x0=0.994, ctrl_y0=0, ctrl_vx0=0, ctrl_vy0=-2.0016, ctrl_T=11.12` | A figure-eight trajectory that loops around the Moon and returns to the vicinity of Earth without propulsion. Used by Apollo missions as an abort trajectory. | Always |
| `btn_reset` | `ui.iconFA('undo')` | Reset all controls to default values | Reset all parameters to default values | Always |

### 3.5. Custom UI Components

| Component ID | Brief description | Role | UI component specification |
|---|---|---|---|
| `comp_lagrange_info` | Display panel showing Lagrange point coordinates and stability. Built with `ui.divV`/`ui.label`; values updated via `textContent`. Shows L1–L5 with x,y coordinates. L4/L5 stability note based on μ vs Routh criterion. Styled via `.cr3bp-app-lagrange-panel`. | Display | Inline |

---

## 4. Result Display Elements

### 4.1. Primary Pipeline Display Elements

| ID | Type | Associated output data | Docking location |
|---|---|---|---|
| `view_orbit` | Datagrok viewer `scatter plot` | `x`, `y` columns, colored by `Velocity` | Right area, top |
| `view_jacobi` | Datagrok viewer `line chart` | `Time`, `Jacobi Constant` columns | Right area, bottom-right |
| `view_zvc` | Datagrok viewer `scatter plot` | ZVC grid DataFrame (`x_grid`, `y_grid`, `Forbidden`) | Right area, bottom-left |
| `view_table` | `DG.TableView` grid | `Time`, `x`, `y`, `vx`, `vy`, `Velocity`, `Jacobi Constant` | Center (default) |
| `comp_lagrange_info` | Custom HTMLElement | `lagrangePoints` from solution | Left panel, below controls |

**Details for `view_orbit`:**

| Property | Value |
|---|---|
| X axis | `x — distance along primary axis` |
| Y axis | `y — distance perpendicular to axis` |
| Color | `Velocity` (velocity magnitude) — linear color scheme: cool tones (low) to warm tones (high) |
| Title | Orbit in Rotating Frame (x, y) |
| Markers | Body 1 (blue, large, at (−μ, 0), label "Earth"), Body 2 (gray, large, at (1−μ, 0), label "Moon"). L1–L5 as small labeled markers. Rendered via `onAfterDrawScene` overlay. |

Lagrange point labels in overlay:

| Marker | Label | Description |
|---|---|---|
| L1 | L1 — inner equilibrium | Between the two bodies |
| L2 | L2 — far-side equilibrium | Beyond the smaller body |
| L3 | L3 — opposite equilibrium | Beyond the larger body |
| L4 | L4 — leading triangle point | Leading equilateral triangle vertex |
| L5 | L5 — trailing triangle point | Trailing equilateral triangle vertex |

**Details for `view_jacobi`:**

| Property | Value |
|---|---|
| X axis | Time |
| Y axis | Jacobi Constant |
| Title | Jacobi Constant Drift (integration accuracy) |

**Details for `view_zvc`:**

| Property | Value |
|---|---|
| X axis | `x_grid` |
| Y axis | `y_grid` |
| Color | `Forbidden` — binary: semi-transparent fill for forbidden regions, transparent for accessible |
| Title | Forbidden Regions (Hill surfaces at current C_J) |
| Note | Uses a separate DataFrame with the grid points. Only forbidden-region points are included (where `2·U(x,y) > C_J(0)`). |

### 4.2. Secondary Task Display Elements

N/A — no secondary tasks.

---

## 5. Layout and UI Element Placement

### 5.1. Control Placement

| Area | Content |
|---|---|
| Left panel | `ui.form` with grouped controls + Lagrange info panel |
| Ribbon (Presets group) | `btn_preset_earth_moon`, `btn_preset_sun_jupiter`, `btn_preset_pluto_charon` |
| Ribbon (Orbits group) | `btn_preset_halo`, `btn_preset_freereturn` |
| Ribbon (Actions group) | `btn_reset` |
| Main area | Grid (default) |

**Structure of left panel form:**

```
ui.h2('Initial Conditions & System Parameters')

ui.h3('System')
  ctrl_mu

ui.h3('Position')
  ctrl_x0
  ctrl_y0

ui.h3('Velocity')
  ctrl_vx0
  ctrl_vy0

ui.h3('Integration')
  ctrl_T

ui.h2('Lagrange Points')
  comp_lagrange_info
```

### 5.2. Display Element Placement

| Element ID | Docking area | Position / ratio |
|---|---|---|
| `form` (left panel) | `DG.DOCK_TYPE.LEFT` relative to root | ratio `0.6` |
| `view_orbit` | `DG.DOCK_TYPE.RIGHT` relative to root | ratio `0.7` |
| `view_zvc` | `DG.DOCK_TYPE.DOWN` relative to orbit node | ratio `0.5` |
| `view_jacobi` | `DG.DOCK_TYPE.RIGHT` relative to ZVC node | ratio `0.5` |
| `view_table` (grid) | Default position (center) | — |

### 5.3. Styles

CSS file: `css/cr3bp.css`. Import: `import '../../css/cr3bp.css'`. All custom classes use the `cr3bp-app-` prefix.

**Static styles:**

| Element | CSS class(es) | Description |
|---|---|---|
| Lagrange info panel | `.cr3bp-app-lagrange-panel` | Background, padding, border-radius, gap |
| Lagrange point entry | `.cr3bp-app-lagrange-entry` | Flex row, monospace values |
| Lagrange stability note | `.cr3bp-app-lagrange-stable` / `.cr3bp-app-lagrange-unstable` | Green (stable) / red (unstable) text |
| Preset button group | `.cr3bp-app-preset-group` | Button group styling within ribbon |

**Dynamic styles:**

| Element | CSS class(es) | Condition | Description |
|---|---|---|---|
| Lagrange stability note | `.cr3bp-app-lagrange-stable` ↔ `.cr3bp-app-lagrange-unstable` | `μ < 0.0385` vs `μ ≥ 0.0385` | L4/L5 stability changes with μ |

---

## 6. User Feedback

### 6.1. Control Tooltips

| Control ID | Tooltip text | Mechanism |
|---|---|---|
| `ctrl_mu` | Ratio of the smaller body's mass to the total. For Earth–Moon, μ ≈ 0.012. Increasing μ enlarges the gravitational sphere of influence of the second body — the topology of forbidden regions changes. | `tooltipText` property |
| `ctrl_x0` | Starting x coordinate of the test body in the rotating frame. Can also be set by clicking on the orbit plot. | `tooltipText` property |
| `ctrl_y0` | Starting y coordinate of the test body in the rotating frame. Can also be set by clicking on the orbit plot. | `tooltipText` property |
| `ctrl_vx0` | Initial velocity in the x direction. Can be set by dragging from the click point on the orbit plot. | `tooltipText` property |
| `ctrl_vy0` | Initial velocity in the y direction. Can be set by dragging from the click point on the orbit plot. | `tooltipText` property |
| `ctrl_T` | Total integration time in dimensionless units. One lunar orbital period ≈ 2π ≈ 6.28 time units. | `tooltipText` property |
| `btn_preset_halo` | A periodic orbit around L1 — used by real space missions (e.g., JWST). Requires precise initial conditions. | `ui.button` tooltip argument |
| `btn_preset_freereturn` | A figure-eight trajectory that loops around the Moon and returns to the vicinity of Earth without propulsion. Used by Apollo missions as an abort trajectory. | `ui.button` tooltip argument |
| `btn_preset_earth_moon` | Set mass parameter to the Earth–Moon system value | `ui.button` tooltip argument |
| `btn_preset_sun_jupiter` | Set mass parameter to the Sun–Jupiter system value | `ui.button` tooltip argument |
| `btn_preset_pluto_charon` | Set mass parameter to the Pluto–Charon system value | `ui.button` tooltip argument |
| `btn_reset` | Reset all parameters to default values | `ui.iconFA` third argument |

### 6.2. Validators as Feedback

| Input ID | Validation source | Description |
|---|---|---|
| `ctrl_mu` | `validate()` from core | Inline hint on invalid value (out of domain) |
| `ctrl_T` | `validate()` from core | Inline hint on invalid value (non-positive) |

Note: `ctrl_x0`, `ctrl_y0`, `ctrl_vx0`, `ctrl_vy0` accept all real numbers — no validation needed beyond the input type constraint.

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
| `val_01` | `mu ≤ 0` | `ctrl_mu` | "Mass parameter μ must be positive" |
| `val_02` | `mu > 0.5` | `ctrl_mu` | "Mass parameter μ must not exceed 0.5 (by convention m₂ ≤ m₁)" |
| `val_03` | `T ≤ 0` | `ctrl_T` | "Integration time must be positive" |

#### Validation Order

```
1. val_01
2. val_02 (checked only if val_01 passed)
3. val_03
```

All rules are independent except val_01/val_02 (both concern `ctrl_mu`).

#### Returned Map Format

```
Map<inputId, errorMessage>
```

Where `inputId = 'ctrl_mu' | 'ctrl_T'`.

### 7.2. Secondary Task Validation

N/A — no secondary tasks.

---

## 8. Main Pipeline

### 8.1. Primary Pipeline

| Step | Description |
|---|---|
| Parameter input | User sets values through 6 controls (or click-to-launch / presets) → `getInputs()` converts to `CR3BPParams` |
| Validation | `validate(inputs)` returns `Map<InputId, string>`. On errors: mark invalid inputs, call `clearResults()` |
| Computation | `solve(inputs)` — synchronous ODE integration via `mrt()`, then Lagrange points, Jacobi constant, ZVC grid |
| Result display | `updateDataFrame(result)` creates new `DG.DataFrame`, assigns to view and viewers; `updateZVCDataFrame(result)` updates the ZVC scatter plot; `updateLagrangePanel(result)` updates `comp_lagrange_info` |

Reactive trigger: `onValueChanged` with debounce 100 ms.

Error behavior on validation failure: clear results (empty DataFrame), show inline validation errors.

Error behavior on computation failure: clear results + `grok.shell.error(msg)`.

### 8.2. Secondary Pipelines

N/A — no secondary tasks.

### 8.3. Common Pipeline Aspects

#### Control Behavior During Computations

| Pipeline / task | Controls blocked | Which controls |
|---|---|---|
| `task_primary` | No | — (computation < 500ms) |

#### Computation Error Handling

| Pipeline / task | Strategy | Notification method |
|---|---|---|
| `task_primary` | Clear results | `grok.shell.error` |

### 8.4. Computation Blocking and Batch Input Updates

| Scenario | Source | Target controls (ID) | Blocked pipelines |
|---|---|---|---|
| Preset: Halo Orbit | `btn_preset_halo` | `ctrl_mu`, `ctrl_x0`, `ctrl_y0`, `ctrl_vx0`, `ctrl_vy0`, `ctrl_T` | Primary (blocked) |
| Preset: Free-Return | `btn_preset_freereturn` | `ctrl_mu`, `ctrl_x0`, `ctrl_y0`, `ctrl_vx0`, `ctrl_vy0`, `ctrl_T` | Primary (blocked) |
| Preset: Earth-Moon | `btn_preset_earth_moon` | `ctrl_mu` | Primary (blocked) |
| Preset: Sun-Jupiter | `btn_preset_sun_jupiter` | `ctrl_mu` | Primary (blocked) |
| Preset: Pluto-Charon | `btn_preset_pluto_charon` | `ctrl_mu` | Primary (blocked) |
| Reset to defaults | `btn_reset` | All `ctrl_*` | Primary (blocked) |
| Click-to-launch | Scatter plot click/drag | `ctrl_x0`, `ctrl_y0`, `ctrl_vx0`, `ctrl_vy0` | Primary (blocked) |

Reactivity mode during batch update:

| Scenario | Reactivity mode |
|---|---|
| All above | `computationsBlocked = true` → writes → `computationsBlocked = false` → single `runPrimary()` |

---

## 9. Reactivity and Dependencies Between Inputs

### 9.1. Dependency Graph

| Source (input ID) | Target (input IDs) | Reaction type | Logic |
|---|---|---|---|
| `ctrl_mu` | `comp_lagrange_info` | Display update | Lagrange points and stability depend on μ |
| Any `ctrl_*` | All viewers | Reactive update | Rerun `task_primary`, update all displays |

No cascading dependencies between controls — all 6 inputs independently trigger the primary pipeline.

### 9.2. Debounce / Throttle

| Input ID | Strategy | Interval (ms) |
|---|---|---|
| All `ctrl_*` | debounce | 100 |

---

## 10. Data Lifecycle

### 10.1. Data Input

Primary method: manual input via controls, click-to-launch on orbit plot, and preset buttons.

Initial state: all controls initialized with defaults (free-return trajectory for Earth-Moon). `task_primary` runs automatically on initialization.

### 10.2. Loading from Resources

N/A — no external data loading.

### 10.3. Results Table Lifecycle

Update strategy: **DataFrame replacement** — on each recomputation, new DataFrames are created and assigned.

Two DataFrames are managed:
- **Trajectory DataFrame**: columns `[Time, x, y, vx, vy, Velocity, Jacobi Constant]` — assigned to `view.dataFrame` and viewers `view_orbit`, `view_jacobi`, `view_table`.
- **ZVC DataFrame**: columns `[x_grid, y_grid, Forbidden]` — assigned to `view_zvc` only.

```
1. Application initialization
   → solve with DEFAULTS → new trajectory DataFrame + ZVC DataFrame
   → DataFrames assigned to views/viewers

2. task_primary completion
   → new trajectory DataFrame from solution arrays
   → new ZVC DataFrame from zvcGrid
   → assigned to respective views/viewers
   → Lagrange point overlay and body markers redrawn

3. Preset button press
   → batch update controls (computationsBlocked)
   → task_primary runs once → new DataFrames

4. btn_reset press
   → batch update controls to defaults
   → task_primary runs → new DataFrames

5. Click-to-launch
   → batch update x0, y0, vx0, vy0
   → task_primary runs → new DataFrames
```

---

## 11. Error Handling Beyond Computations

| Error type | Strategy | Notification method |
|---|---|---|
| ODE solver error (e.g., singularity at r₁=0 or r₂=0) | Clear results, show message | `grok.shell.error` |
| ODE solver stiffness / non-convergence | Clear results, suggest reducing T or changing initial conditions | `grok.shell.warning` |
| Newton's method failure for L1/L2/L3 | Use fallback approximate values from Hill's approximation | Silent fallback — no user notification |

---

## 12. Subscriptions and Resource Management

### 12.1. Event Subscriptions

| Subscription | Event | Cleanup mechanism |
|---|---|---|
| View removal | `grok.events.onViewRemoved` | `sub.unsubscribe()` in cleanup handler |
| Orbit plot click | `view_orbit` canvas `mousedown`/`mousemove`/`mouseup` | Removed in cleanup handler |
| Orbit plot overlay | `view_orbit.onAfterDrawScene` | `sub.unsubscribe()` in cleanup handler |
| Control changes | `ctrl_*.onValueChanged` | `sub.unsubscribe()` in cleanup handler |

All subscriptions collected in `subs[]` and unsubscribed when the view is removed.

### 12.2. Worker Termination

N/A — no web workers.

---

## 13. Application Closure

On view close (`grok.events.onViewRemoved`), the coordinator performs:

- [x] All event subscriptions unsubscribed (`subs[]` iteration)
- [x] Pending debounce timer cleared (`clearTimeout(debounceTimer)`)
- [x] No workers to terminate
- [x] No open dialogs to close

Closure handler: `grok.events.onViewRemoved.subscribe(...)`.

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
| `src/tests/cr3bp-math-tests.ts` | Math: CR3BP RHS, Math: Lagrange Points, Math: Jacobi Constant, Math: ZVC Grid, Math: MRT solver | ~20 | ODE RHS verification, Lagrange point computation, Jacobi constant conservation, ZVC grid correctness, MRT solver reference problems |
| `src/tests/cr3bp-validation-tests.ts` | API: Validation | ~10 | Validation rules for μ and T, defaults, boundaries, multiple errors |

Tests are run via `grok test` (entry point: `src/package-test.ts`).

### 15.2. Inputs

| Category | Coverage | Description |
|---|---|---|
| Boundary values | `v_bnd_mu_low` (`mu=0.001`), `v_bnd_mu_high` (`mu=0.5`), `v_bnd_T` (`T=1`) | Allowed boundary values |
| Invalid values | `v_01`–`v_03` (6 tests) | Zero, negative, and out-of-range for μ and T |
| Valid defaults | `v_def` | All defaults pass validation |
| Multiple errors | `v_multi` | `mu=0, T=0` → 2 errors |

### 15.3. Mathematical Verification

#### Level 1 Verification

**Formula/equation verification:**

| Test category | Test count | What is verified | Reference source |
|---|---|---|---|
| Math: CR3BP RHS | 5 | ODE RHS: concrete `(μ, x, y, vx, vy)` → expected `(dx/dt, dy/dt, dvx/dt, dvy/dt)` | Manual calculation |
| Math: Lagrange Points | 4 | L4/L5 analytical values; L1/L2/L3 numerical values for Earth-Moon μ | Analytical formulas; Szebehely (1967) |
| Math: Jacobi Constant | 2 | `C_J` computation for known state vectors | Manual calculation |

**Output property verification:**

| Test category | Test count | Properties verified |
|---|---|---|
| Math: Jacobi conservation | 2 | `|C_J(t) − C_J(0)| < 1e-6` for all t, checked at default params and a second parameter set |
| Math: Lagrange equilibrium | 3 | L1/L2/L3: `|f(x_L)| < 1e-10` (equilibrium condition); L4/L5: accelerations ≈ 0 at triangular points |
| Math: ZVC consistency | 1 | All grid points marked forbidden satisfy `2·U(x,y) > C_J(0)`, all accessible satisfy `2·U(x,y) ≤ C_J(0)` |

#### Level 2 Verification

**Numerical method verification:**

| Test category | Test count | Reference problems | Tolerance | Source |
|---|---|---|---|---|
| Math: MRT solver | 2 | Known periodic orbit (Halo near L1 with known period), Two-body limit (μ→0, compare to Keplerian ellipse) | Endpoint error < 0.01 | Richardson (1980); Keplerian analytics |

**Convergence verification:**

| Status | Description |
|---|---|
| Implemented | Solve Halo orbit preset with tolerance 1e-6 and 1e-8, verify endpoint discrepancy decreases |

**Asymptotic/equilibrium behavior:**

| Status | Description |
|---|---|
| Implemented | Start at L4 with μ=0.01 (stable), integrate for T=100, verify body remains near L4 (`|Δx|`, `|Δy|` < 0.01) |

#### Validation Tests

| Test ID | Rule | Input data | Expected result |
|---|---|---|---|
| `v_01a` | val_01 | `mu=0` | Error on `ctrl_mu` |
| `v_01b` | val_01 | `mu=-0.1` | Error on `ctrl_mu` |
| `v_02a` | val_02 | `mu=0.6` | Error on `ctrl_mu` |
| `v_02b` | val_02 | `mu=1.0` | Error on `ctrl_mu` |
| `v_03a` | val_03 | `T=0` | Error on `ctrl_T` |
| `v_03b` | val_03 | `T=-10` | Error on `ctrl_T` |
| `v_def` | all defaults | DEFAULTS | `errors.size = 0` |
| `v_bnd_mu_low` | valid boundary | `mu=0.001` | No error on `ctrl_mu` |
| `v_bnd_mu_high` | valid boundary | `mu=0.5` | No error on `ctrl_mu` |
| `v_multi` | multiple | `mu=0, T=0` | 2 errors |

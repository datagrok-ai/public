# CLAUDE.md

## Overview

**Diff Studio** is a Datagrok package that provides in-browser tools for solving initial value problems (IVP) for systems of ordinary differential equations (ODEs). It implements multiple numerical methods for solving both stiff and non-stiff ODEs directly in the browser, including the automatic stiffness-detecting LSODA method, Rosenbrock-Wanner implicit methods, Runge-Kutta explicit methods, and Adams-Bashforth multistep methods.

The package is accessible via **Apps > Compute > Diff Studio** in the Datagrok platform.

## Key Dependencies

- **diff-grok** (v1.2.0+) - Core ODE solver library implementing CVODE, LSODA, MRT, ROS3PRw, ROS34PRw, RK3, RK4, RKDP, AB4, and AB5 methods
- **@datagrok-libraries/compute-utils** - Provides sensitivity analysis and fitting views
- **CodeMirror 6** - Code editor for IVP formula editing
- **@datagrok-libraries/test** - Testing utilities

## Architecture

### Core Components

**app.ts (DiffStudio class)**
- Main UI application for the solver
- Handles IVP code editing with CodeMirror
- Manages templates, library examples, and user models
- Integrates with sensitivity analysis and parameter fitting
- File preview and browser integration
- ~2300+ lines - the heart of the user interface

**solver-tools.ts**
- Thin wrapper around `diff-grok` methods (lsoda, mrt, ros3prw, ros34prw, rk3, rk4, rkdp, ab4, ab5)
- `solveDefault()` - uses ROS34PRw by default
- `solveIVP()` - customizable solver with options
- `getMethod()` - resolves method name from solver options
- Converts solver output to Datagrok DataFrames

**scripting-tools.ts**
- Parser for IVP declarative format (`.ivp` files)
- Generates JavaScript script code from IVP specifications
- `getIVP()` - parses IVP text to structured format
- `getScriptLines()` - converts IVP to executable JS code
- Handles all control blocks (#name, #equations, #argument, #inits, etc.)

**model.ts**
- Represents a complete differential equations model
- Manages model execution with Diff Studio UI
- Handles dock layout ratios for inputs/graphs

**utils.ts**
- `error()` - Calculates max absolute deviation between DataFrames
- `unusedFileName()` - Generates unused IVP file names

**callbacks/** directory
- `callback-base.ts` - Base callback interface
- `callback-tools.ts` - Callback manager, creates callbacks based on solver options
- `iter-checker-callback.ts` - Iteration limit checking
- `time-checker-callback.ts` - Computation time limit checking

**export/** directory
- `export-dialog.ts` - Interactive LaTeX/Markdown export dialog for the current model
- Uses `convertIvpToLatex` from `diff-grok` to render the IVP text
- Live-preview is a read-only CodeMirror 6 editor; language (`stex` for LaTeX via `StreamLanguage`, `markdown()` for Markdown) is swapped in-place through a `Compartment` when the format input changes
- Offers clipboard copy and file download; LaTeX output can be wrapped into a standalone compilable document

### Templates and Examples

**templates.ts** - Three starter templates:
- Basic - simplest model
- Advanced - includes expressions, constants, parameters, tolerance
- Extended - demonstrates annotation features for UI generation

**use-cases.ts** - Built-in examples:
- Chemical reactions (mass-action kinetics)
- Robertson's model (stiff equations)
- Fermentation kinetics
- PK (pharmacokinetic)
- PK-PD (pharmacokinetic-pharmacodynamic)
- Acid production
- Nimotuzumab (population pharmacokinetic)
- Bioreactor processes
- Pollution (air pollution model, 25 reactions)

**demo/** directory - Standalone model demonstrations:
- `acid-production.ts` - Gluconic acid production model
- `ball-flight.ts` - Ball trajectory simulation
- `bioreactor.ts` - Controlled fab-arm exchange mechanism
- `pk-pd.ts` - PK-PD simulation demo
- `pollution.ts` - Air pollution model

### Files Structure

```
files/
  templates/          # .ivp template files
    basic.ivp
    advanced.ivp
    extended.ivp
  library/            # Example .ivp files
    chem-react.ivp
    robertson.ivp
    fermentation.ivp
    pk.ivp
    pk-pd.ivp
    ga-production.ivp
    nimotuzumab.ivp
    bioreactor.ivp
    pollution.ivp
    energy-n-control.ivp
  icons/              # Model icons
  ball-flight-trajectory.csv
```

## IVP File Format

The `.ivp` format is a declarative specification for initial value problems. Control blocks:

```python
#name: Model Name
#description: Model description
#equations:
  dx/dt = expression
  dy/dt = expression
#argument: t
  initial = 0
  final = 10
  step = 0.01
#inits:
  x = 1.0
  y = 2.0
#constants:        # optional
  k = 3.14
#parameters:       # optional
  P = 1 {min: 0; max: 5}
#expressions:      # optional
  aux = x + y
#output:           # optional
  Time = t
  X = x {caption: Position}
#meta.solver: {method: 'ros34prw'; maxIterations: 10000; maxTimeMs: 5000}
#meta.inputs: table {choices: OpenFile("...")}
```

Annotations (in `{...}`) control UI generation: `caption`, `category`, `min`, `max`, `step`, `units`, etc.

## Package Functions

**Main app function:**
- `runDiffStudio()` - Launches Diff Studio application

**File handlers:**
- `ivpFileHandler()` - Handles .ivp file double-click
- `previewIvp()` - File previewer for .ivp files

**Solver functions (programmatic API):**
- `solve(problem: ODEs)` - Solve using default options
- `solveEquations(problem: ODEs, options: SolverOptions)` - Solve with custom options
- `solveODE(problem: string)` - Parse and solve IVP from string

**Model functions:**
- `acidProduction()` - Gluconic acid production model
- `ballFlight()` - Ball trajectory model
- `Bioreactor()` - Bioreactor model
- `pkPdNew()` - PK-PD model
- `pollution()` - Air pollution model
- `runModel()` - Run model with Diff Studio UI

**Utility functions:**
- `serializeEquations(problem: string)` - Parse IVP to object
- `odesToCode(serialization: IVP)` - Generate JS code from IVP

**Demos:**
- `runDiffStudioDemo()` - Interactive demo
- `demoSimPKPD()` - PK-PD demo
- `demoBioreactor()` - Bioreactor demo

## Testing

Tests are organized in `src/tests/`:

**numerical-methods-tests.ts**
- Correctness tests - validate solver accuracy against known solutions
- Performance tests - benchmark solver speed
- Uses `corrProbs` and `perfProbs` from `diff-grok`

**features-tests.ts**
- Tests for IVP format features

**platform-funcs-tests.ts**
- Tests for platform function integration

**pipeline-tests.ts**
- End-to-end pipeline tests

**demo-models-tests.ts**
- Tests for demo model functions

**parser-tests.ts**
- Tests for IVP parser

**test-utils.ts**
- Shared test utilities

## Numerical Methods

Diff Studio implements the following methods from `diff-grok`:

**Automatic stiffness-detecting methods:**
- **CVODE** - `method: 'cvode'` - variable-order, variable-step BDF solver from SUNDIALS v7.5.0 port (Hindmarsh et al., 2005); uses dense direct linear solver (LU decomposition) with warmup strategy for extremely stiff problems
- **LSODA** - `method: 'lsoda'` - variable-order Nordsieck-based solver with automatic switching between Adams (non-stiff) and BDF (stiff)

**Implicit methods (for stiff ODEs) - Rosenbrock-Wanner type:**
- **MRT** (Modified Rosenbrock Triple) - `method: 'mrt'`
- **ROS3PRw** - `method: 'ros3prw'`
- **ROS34PRw** (default) - `method: 'ros34prw'`

**Explicit methods (for non-stiff ODEs) - Runge-Kutta type:**
- **RK3** (Bogacki-Shampine 3(2)) - `method: 'rk3'`
- **RK4** (Runge-Kutta-Fehlberg 4(5)) - `method: 'rk4'`
- **RKDP** (Dormand-Prince 5(4)) - `method: 'rkdp'`

**Explicit methods (for non-stiff ODEs) - Adams-Bashforth type:**
- **AB4** (predictor-corrector of order 4) - `method: 'ab4'`
- **AB5** (predictor-corrector of order 5) - `method: 'ab5'`

ROS34PRw is the default method. CVODE and LSODA are recommended as general-purpose solvers that auto-detect stiffness.

## Important Constants

**constants.ts** - Core constants (DF_NAME, CONTROL_EXPR, MAX_LINE_CHART, etc.)
**ui-constants.ts** - UI-specific constants (TITLE, HINT, ERROR_MSG, DOCK_RATIO, etc.)

## Error Handling

**error-utils.ts** provides:
- `ModelError` - Custom error class for model errors
- `showModelErrorHint()` - Display error hints to users
- Error message generators for common issues

## Key Architecture Notes

- The IVP parser in `scripting-tools.ts` is the bridge between declarative `.ivp` format and executable code
- `app.ts` manages complex UI state with CodeMirror, file browsing, model history, and analysis integration
- Web workers are available via `getIvp2WebWorker()` from `diff-grok` for background solving
- The package uses Datagrok's `@grok.decorators` for function registration (not plain JSDoc comments)
- Models can be saved to "Model Hub" for sharing and reuse
- Sensitivity analysis and parameter fitting are integrated via `@datagrok-libraries/compute-utils`
- Tree-click previews reuse a single `DiffStudio` instance via the static `sharedPreview` slot. The first click runs a cold path (`getStatePreview`/`getFilePreview`); subsequent clicks call `applyState`/`applyFile`, which swap the CodeMirror doc via `EditorState.create` and re-solve without remounting. **Any new transient per-model field on `DiffStudio` must also be drained in `resetForReuse()`** — otherwise it leaks between consecutive previews.
- Any path that opens a full-view (creates a `DiffStudio` and calls `grok.shell.addView` on its `solverView`) must call `DiffStudio.releaseSharedPreview()` first, so the static slot does not pin a stale preview after the user navigates away. `runSolverApp` already does this; new entry points should follow the same pattern.

## Styling conventions

- **Put styles in CSS, not inline.** Define classes in `css/app-styles.css` and apply them with `element.classList.add('...')`. Avoid `element.style.cssText = '...'` and `element.style.prop = '...'` in TypeScript, except for values that must be computed at runtime (e.g. toggling `display: none`, setting a width derived from a measurement).
- Class names follow the `diff-studio-*` prefix for namespacing.
- When adding a new UI element, first check `css/app-styles.css` for a reusable class; only introduce a new one if none fits.

## Code style (functions)

- **Every function gets a JSDoc comment.** Place a `/** ... */` block immediately above each function (top-level functions, exported functions, inner arrow functions, and methods) describing what it does. Keep it to one or two sentences focused on the WHY / contract, not a line-by-line restatement of the body.
- **Mark the closing brace of long functions.** When a function body exceeds 20 lines, append `// <functionName>` to its closing brace so readers can see which function is ending. Example: `} // showExportDialog`. Do not add this marker to shorter functions.

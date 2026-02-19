# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**Diff Studio** is a Datagrok package that provides in-browser tools for solving initial value problems (IVP) for systems of ordinary differential equations (ODEs). It implements Rosenbrock-Wanner numerical methods for solving both stiff and non-stiff ODEs directly in the browser.

The package is accessible via **Apps > Compute > Diff Studio** in the Datagrok platform.

## Key Dependencies

- **diff-grok** (v1.0.8+) - Core ODE solver library implementing MRT, ROS3PRw, and ROS34PRw methods
- **@datagrok-libraries/compute-utils** - Provides sensitivity analysis and fitting views
- **CodeMirror 6** - Code editor for IVP formula editing
- **@datagrok-libraries/test** - Testing utilities

## Build Commands

```bash
# Standard build pipeline
npm run build              # grok api && grok check --soft && webpack

# Development/publishing
npm run debug-odes         # webpack && grok publish (default server)
npm run debug-odes-dev     # webpack && grok publish dev
npm run debug-odes-local   # webpack && grok publish local
npm run release-odes       # webpack && grok publish --release

# Testing
npm run test               # grok test (against default server)

# Individual steps
npm run build-odes         # webpack only
grok api                   # Generate package.g.ts and package-api.ts
grok check --soft          # Validate package
```

## Architecture

### Core Components

**app.ts (DiffStudio class)**
- Main UI application for the solver
- Handles IVP code editing with CodeMirror
- Manages templates, library examples, and user models
- Integrates with sensitivity analysis and parameter fitting
- File preview and browser integration
- ~2500+ lines - the heart of the user interface

**solver-tools.ts**
- Thin wrapper around `diff-grok` methods (mrt, ros3prw, ros34prw)
- `solveDefault()` - uses ROS34PRw by default
- `solveIVP()` - customizable solver with options
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

**callbacks/** directory
- `callback-base.ts` - Base callback interface
- `iter-checker-callback.ts` - Iteration limit checking
- `time-checker-callback.ts` - Computation time limit checking

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
- `ball-flight.ts` - Ball trajectory simulation
- `pk-pd.ts` - PK-PD simulation demo
- `bioreactor.ts` - Controlled fab-arm exchange mechanism

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
- `ballFlight()` - Ball trajectory model
- `pkPdNew()` - PK-PD model
- `Bioreactor()` - Bioreactor model
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

Run specific test categories:
```bash
grok test --category "Correctness"
grok test --category "Performance"
grok test --test "MRT"
grok test --gui              # Visual debugging
```

## Development Workflow

1. **Edit IVP format parser** - modify `scripting-tools.ts`
2. **Change UI behavior** - modify `app.ts` (DiffStudio class)
3. **Adjust solver integration** - modify `solver-tools.ts`
4. **Add templates/examples** - add to `templates.ts`, `use-cases.ts`, or `files/`
5. **Update model demos** - modify files in `demo/`

After changes:
```bash
npm run build        # Build and validate
npm run debug-odes   # Publish to test server
```

## Numerical Methods

Diff Studio implements three Rosenbrock-Wanner methods from `diff-grok`:

- **MRT** (Modified Rosenbrock Triple) - `method: 'mrt'`
- **ROS3PRw** - `method: 'ros3prw'`
- **ROS34PRw** (default) - `method: 'ros34prw'`

All methods handle stiff and non-stiff equations. ROS34PRw is the default for best accuracy/performance balance.

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

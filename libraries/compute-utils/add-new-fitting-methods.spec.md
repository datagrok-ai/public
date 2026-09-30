# Spec: New optimization methods in the Fitting View (compute-utils ← sci-comp)

> Development specification. Nothing is implemented or committed at this stage.

## 1. Context

The Fitting View in the `compute-utils` library (single-objective model-parameter fitting) currently
offers a single optimization method — **Nelder-Mead** (a self-contained implementation inside
compute-utils). Meanwhile the `@datagrok-libraries/sci-comp` library (same author, published on npm)
already provides a ready-made, test-covered set of optimizers behind a single interface. The goal is
to **reuse** the sci-comp optimizers and surface them in the Fitting View without reimplementing the
algorithms.

The problem to solve: the two libraries have different optimizer contracts (function signatures,
settings format, result format, bounds handling), so we need an adapter and an extensible method
registration layer.

Outcome: the Fitting View gains a method selector (`L-BFGS-B`, `PSO`, `L-BFGS`, `Adam` in addition to
`Nelder-Mead`), each method with its own settings in the UI, with correct handling of parameter bounds.

## 2. Decisions made (from discussion)

| Question | Decision |
|---|---|
| Which methods to add | **L-BFGS-B, PSO, L-BFGS, Adam** (Nelder-Mead stays its own). Gradient Descent is NOT included — its unfinished scaffolding is removed. |
| Parameter bounds | **Native + smooth penalty**: L-BFGS-B — native box bounds; L-BFGS/Adam/PSO — sci-comp smooth quadratic penalty (`boxConstraints`). Requires passing bounds into the optimizer. |
| Architecture | **Method registry** (declarative `METHOD → descriptor`), modeled on `OptimizeManager` from the multi-objective `optimization-view`. |
| Execution (standard path) | **Main-thread only** for the new methods. The standard path's worker pool stays NM-only (it uses sync codegen); worker support there is a separate follow-up. |
| Diff Studio (IVP models) | **In scope.** The new methods work for Diff Studio too — directly inside its worker (`diff-studio/workers/basic.ts` → `fit()`), because it uses an **async** optimizer → codegen / sync twins are NOT needed. |
| Iteration budgets | **Conservative** defaults in the descriptors (reference — the current NM = 50 iterations), the user raises them manually. NOT the sci-comp defaults (1000/10000). |
| Penalty coefficient μ | **Expose μ as a method setting** (for L-BFGS/Adam/PSO), default 1000. |
| Gradient cost | **Cap with a default** `maxFunctionEvaluations`/`maxIterations` so gradient methods don't hang on an expensive model. |
| Nelder-Mead | **Keep the self-contained** compute-utils NM (preserves the standard path's worker fast-path + codegen). New methods go through sci-comp. |

## 3. How it works today (brief, with anchors)

Fitting engine: [function-views/src/fitting/](function-views/src/fitting/)

- **Optimizer contract** `IOptimizer` — [optimizer-misc.ts:34-40](function-views/src/fitting/optimizer-misc.ts#L34-L40):
  ```ts
  (objectiveFunc: (x: Float64Array) => Promise<number|undefined>,
   paramsInitial: Float64Array,
   settings: Map<string, number>,
   threshold?: number) => Promise<Extremum>
  ```
  Returns `Extremum = {point, cost, iterCosts, iterCount}` ([optimizer-misc.ts:6-11](function-views/src/fitting/optimizer-misc.ts#L6-L11)).
  Bounds are NOT passed to the optimizer: `objectiveFunc` returns `undefined` outside the bounds, and NM
  substitutes the penalty `costOutside = 2×(initial cost)` ([optimizer-nelder-mead.ts:137](function-views/src/fitting/optimizer-nelder-mead.ts#L137)).
- **Nelder-Mead** — `optimizeNM: IOptimizer` + a declarative settings descriptor
  `nelderMeadSettingsOpts: Map<string, Setting>` ([optimizer-nelder-mead.ts:7-64, 121](function-views/src/fitting/optimizer-nelder-mead.ts#L7-L64)).
  Type `Setting = {default, min, max, caption, tooltipText, inputType}` ([optimizer-misc.ts:13-20](function-views/src/fitting/optimizer-misc.ts#L13-L20)).
- **Dispatcher** `performNelderMeadOptimization` — [optimizer.ts:9-61](function-views/src/fitting/optimizer.ts#L9-L61): picks worker/main-arm; **hard-wired to NM**.
- **Main-arm** `MainExecutor.run` — [worker/executor.ts:70-111](function-views/src/fitting/worker/executor.ts#L70-L111): samples starting points (`sampleSeeds`), loops calling `optimizeNM(...)` directly ([executor.ts:90](function-views/src/fitting/worker/executor.ts#L90)). `canHandle` ([executor.ts:57-66](function-views/src/fitting/worker/executor.ts#L57-L66)) decides about the worker, knows nothing about the method.
- **Public API** `optimizer-api.ts`: `OptimizerParams`, `resolveParams` (merges only `nelderMeadSettingsOpts`), `runRawOptimizer` (always calls NM), `runOptimizerFinalized` — [optimizer-api.ts:14-113](function-views/src/fitting/optimizer-api.ts#L14-L113).
- **Bounds**: `getAccData(inputs)` (in `bounds-checker.ts`) splits the config into `constValues / nonFormulaBounds / formulaBounds`; the varied-parameter order is `[...nonFormulaBounds, ...formulaBounds]` (see [optimizer-sampler.ts:45-82](function-views/src/fitting/optimizer-sampler.ts#L45-L82)). Each entry carries `boundsIdx` — its position in the parameter vector.
- **UI** (single-objective) — [function-views/src/fitting-view.ts](function-views/src/fitting-view.ts):
  - method dropdown, `items: [METHOD.NELDER_MEAD]` — [fitting-view.ts:689-694](function-views/src/fitting-view.ts#L689-L694);
  - auto-generation of settings inputs from the NM descriptor — [fitting-view.ts:902-937](function-views/src/fitting-view.ts#L902-L937);
  - validation against the descriptor — [fitting-view.ts:940-949](function-views/src/fitting-view.ts#L940-L949);
  - hard guard `method !== NELDER_MEAD → throw` at run time — [fitting-view.ts:1621-1622](function-views/src/fitting-view.ts#L1621-L1622);
  - Gradient Descent scaffolding (`gradDescentSettings`, `generateGradDescentSettingsInputs`, a validator that `return false`) — dead code, to be removed.
- **Reference extensible pattern** (multi-objective): `optManagers = Map<METHOD, OptimizeManager>` and the abstract `OptimizeManager` (`getInputs / areSettingsValid / perform / visualize`) — [function-views/src/optimization-view.ts:366-368](function-views/src/optimization-view.ts#L366-L368) and `multi-objective-optimization/optimize-manager.ts`.
- **METHOD enum** — [constants.ts:4-12](function-views/src/fitting/constants.ts#L4-L12).

## 4. What sci-comp provides

[../sci-comp/src/optimization/single-objective/](../sci-comp/src/optimization/single-objective/), npm `@datagrok-libraries/sci-comp@0.11.0`.

- **Base class** `Optimizer<S extends CommonSettings>` — [optimizer.ts:23](../sci-comp/src/optimization/single-objective/optimizer.ts#L23). Public: `minimize/maximize` (sync) and `minimizeAsync/maximizeAsync`. Wraps `constraints` into a penalty itself ([optimizer.ts:140-151](../sci-comp/src/optimization/single-objective/optimizer.ts#L140-L151)).
- **Types** — [types.ts](../sci-comp/src/optimization/single-objective/types.ts): `AsyncObjectiveFunction = (x)=>Promise<number>`; `OptimizationResult = {point, value, iterations, converged, costHistory}`; `CommonSettings = {maxIterations?, tolerance?, onIteration?, constraints?, penaltyOptions?}`; `IterationCallback = (state)=>boolean|void` (return `true` to stop); `Constraint = {type:'ineq'|'eq', fn}`.
- **Penalty layer** `boxConstraints(lower, upper): Constraint[]` — converts a box into inequalities (smooth quadratic penalty; ±∞ entries are harmless). `applyPenalty` — μ defaults to 1000.
- **Registry** `getOptimizer(name)`, `listOptimizers()`, `registerOptimizer`.
- **Classes**: `NelderMead`, `PSO`, `GradientDescent`, `Adam`, `LBFGS`, `LBFGSB` (all under the `singleObjective` namespace from the package root) — [index.ts](../sci-comp/src/optimization/single-objective/index.ts).
- **L-BFGS-B** `LBFGSBSettings`: `bounds?: {lower?, upper?}` (native box bounds, scalar/array/±∞), `gradFn?` (analytic gradient, optional), `gradTolerance` (pgtol), `maxFunctionEvaluations`, `lineSearch{ftol,gtol,xtol,maxSteps}`, `historySize`.
- **PSO** `PSOSettings`: `swarmSize, inertia, cognitive, social, searchRange?:{lower,upper}, seed?, noImprovementLimit`.
- **L-BFGS** `LBFGSSettings`: `historySize, gradTolerance, c1, maxLineSearchSteps, initialStepSize, finiteDiffStep`.
- **Adam** `AdamSettings`: `learningRate, beta1, beta2, epsilon, finiteDiffStep, gradTolerance, maxGradNorm, noImprovementLimit`.

Every single-objective optimizer has an async variant — matching exactly what `IOptimizer` needs.

## 5. Target architecture

### 5.1 Dependency
Add to [package.json](package.json) under `dependencies`:
`"@datagrok-libraries/sci-comp": "workspace:^"`. Both libraries are in the same pnpm workspace,
Turborepo builds sci-comp first. Import: `import {singleObjective} from '@datagrok-libraries/sci-comp';`
(`singleObjective.LBFGSB`, `singleObjective.boxConstraints`, `singleObjective.Optimizer`, settings types).

### 5.2 Adapter `Optimizer` (sci-comp) → `IOptimizer` (compute-utils)
New file `fitting/optimizer-sci-comp.ts`. Extend the contract with bounds (optional 5th parameter,
NM ignores it — backward compatible):

```ts
// optimizer-misc.ts — interface extension
export type OptimizerBounds = { lower: Float64Array; upper: Float64Array }; // ±Infinity where unbounded
export interface IOptimizer {
  (objectiveFunc, paramsInitial, settings, threshold?, bounds?: OptimizerBounds): Promise<Extremum>;
}
```

Generic adapter:
```ts
export function adaptOptimizer<S extends CommonSettings>(
  makeOpt: () => singleObjective.Optimizer<S>,
  buildSettings: (m: Map<string, number>, bounds?: OptimizerBounds) => S,
): IOptimizer {
  return async (objectiveFunc, x0, settings, threshold, bounds) => {
    const opt = makeOpt();
    const s = buildSettings(settings, bounds);
    // undefined (out of bounds) → penalty, mirroring NM's costOutside
    const penalty = 2 * ((await objectiveFunc(x0)) ?? Infinity);
    const fn = async (x: Float64Array) => (await objectiveFunc(x)) ?? penalty;
    if (threshold != null) {
      const prev = s.onIteration;
      s.onIteration = (st) => (st.bestValue <= threshold) || (prev?.(st) === true);
    }
    const r = await opt.minimizeAsync(fn, x0, s);
    return {point: r.point, cost: r.value, iterCosts: Array.from(r.costHistory), iterCount: r.iterations};
  };
}
```
Result mapping: `value→cost`, `costHistory→iterCosts` (for the "Fitting profile" chart), `iterations→iterCount`.

### 5.3 Per-method settings descriptors
In the same file — one `Map<string, Setting>` per method (`psoSettingsOpts`, `adamSettingsOpts`,
`lbfgsSettingsOpts`, `lbfgsbSettingsOpts`), analogous to `nelderMeadSettingsOpts`: each key is a
sci-comp settings field with `caption/tooltipText/min/max/inputType/default`. Plus, per method, a
`buildSettings(m, bounds)` function that translates the flat `Map<string,number>` into typed `*Settings`
and injects the bounds (see 5.4). We limit ourselves to a practically useful subset of fields
(e.g. for L-BFGS-B: `maxIterations`, `gradTolerance`, `historySize`, `maxFunctionEvaluations`;
for PSO: `swarmSize`, `inertia`, `cognitive`, `social`, `maxIterations`; etc.).

Per the decisions made:
- **Conservative iteration defaults** (not sci-comp's 1000/10000): moderate per-method values
  (reference — the current NM `maxIter`=50), with the ability to raise them in the UI.
- **`maxFunctionEvaluations` is capped by default** for gradient methods (L-BFGS/Adam/GD) and L-BFGS-B,
  since the numerical gradient ≈ `2·dim` model runs per iteration.
- **μ (penalty coefficient) is a separate settings field** for L-BFGS/Adam/PSO (default 1000), mapped
  to `penaltyOptions.mu` (see 5.4). For L-BFGS-B with native bounds, μ is not applied (no penalty needed).

### 5.4 Passing bounds (native + smooth penalty)
New helper (next to the sampler, reuses `getAccData`): builds `OptimizerBounds` for the varied
parameters in the same order as the parameter vector:
- `nonFormulaBounds` → finite `lower/upper` at position `boundsIdx`;
- `formulaBounds` → `-Infinity/+Infinity` (natively unbounded; for these the protective undefined
  cutoff from `objectiveFunc` remains).

Application in `buildSettings`:
- **L-BFGS-B**: `s.bounds = {lower, upper}` — native geometric projection.
- **L-BFGS / Adam / PSO**: `s.constraints = singleObjective.boxConstraints(lower, upper)` — smooth
  quadratic penalty (the numerical gradient stays meaningful near the boundary); `s.penaltyOptions = {mu}`
  — μ from the method settings. For PSO additionally `s.searchRange = {lower, upper}` (particle scatter
  within the box; finite bounds, otherwise fall back).

### 5.5 Method registry
New file `fitting/optimizer-registry.ts`:
```ts
export type OptimizerDescriptor = {
  method: METHOD;
  optimizer: IOptimizer;
  settingsOpts: Map<string, Setting>;
  supportsWorker: boolean;   // NM: true, others: false
  wantsBounds: boolean;      // whether to pass OptimizerBounds
};
export const OPTIMIZERS = new Map<METHOD, OptimizerDescriptor>([
  [METHOD.NELDER_MEAD, {optimizer: optimizeNM, settingsOpts: nelderMeadSettingsOpts, supportsWorker: true,  wantsBounds: false, ...}],
  [METHOD.LBFGSB,      {optimizer: adaptOptimizer(() => new singleObjective.LBFGSB(), buildLbfgsbSettings), settingsOpts: lbfgsbSettingsOpts, supportsWorker: false, wantsBounds: true,  ...}],
  [METHOD.PSO,         {optimizer: adaptOptimizer(() => new singleObjective.PSO(),    buildPsoSettings),    settingsOpts: psoSettingsOpts,    supportsWorker: false, wantsBounds: true,  ...}],
  [METHOD.LBFGS,       {optimizer: adaptOptimizer(() => new singleObjective.LBFGS(),  buildLbfgsSettings),  settingsOpts: lbfgsSettingsOpts,  supportsWorker: false, wantsBounds: true,  ...}],
  [METHOD.ADAM,        {optimizer: adaptOptimizer(() => new singleObjective.Adam(),   buildAdamSettings),   settingsOpts: adamSettingsOpts,   supportsWorker: false, wantsBounds: true,  ...}],
]);
```
In [constants.ts](function-views/src/fitting/constants.ts): add to `METHOD`
`LBFGSB='L-BFGS-B'`, `PSO='PSO'`, `LBFGS='L-BFGS'`, `ADAM='Adam'`; **remove** `GRAD_DESC`; update `methodTooltip`.

### 5.6 Engine changes
- `optimizer.ts`: rename/generalize `performNelderMeadOptimization` → `performOptimization`, accept the
  selected `OptimizerDescriptor` (or `method`); pass `optimizer`/`wantsBounds`/`bounds` into `ExecutorArgs`.
- `worker/executor.ts`:
  - `ExecutorArgs` += `optimizer: IOptimizer`, `bounds?: OptimizerBounds` (or `method`).
  - `MainExecutor.run` — when `wantsBounds`, compute `OptimizerBounds` (helper from 5.4) and call
    `args.optimizer(objectiveFunc, params[i], settings, threshold, bounds)` instead of `optimizeNM` ([executor.ts:90](function-views/src/fitting/worker/executor.ts#L90)).
  - `canHandle` — return `false` if the method is not `supportsWorker` (gate on NM), so the worker arm never
    receives a method it doesn't know.

### 5.7 API changes (`optimizer-api.ts`)
- `OptimizerParams` += `method?: METHOD` (defaults to `NELDER_MEAD`).
- `resolveParams` — merge defaults from the selected method's `settingsOpts` (via the registry) rather than
  hard-coded `nelderMeadSettingsOpts` ([optimizer-api.ts:52-55](function-views/src/fitting/optimizer-api.ts#L52-L55)).
- `runRawOptimizer` — take the optimizer/flags from the registry and call `performOptimization`.

### 5.8 UI changes (`fitting-view.ts`)
- `methodInput.items` = `[...OPTIMIZERS.keys()]` ([fitting-view.ts:690](function-views/src/fitting-view.ts#L690)).
- Replace `generateNelderMeadSettingsInputs`/`generateGradDescentSettingsInputs` with a single
  `generateSettingsInputs(method)` that iterates the descriptor's `settingsOpts` (generalizing
  [fitting-view.ts:902-937](function-views/src/fitting-view.ts#L902-L937)).
- Store settings values as `Map<METHOD, Map<string, number>>` (generalizing the single `nelderMeadSettings`).
- `areMethodSettingsCorrect` — validate against the selected method's `settingsOpts` (generalizing [fitting-view.ts:940-949](function-views/src/fitting-view.ts#L940-L949)).
- Remove the `method !== NELDER_MEAD` guard ([fitting-view.ts:1621-1622](function-views/src/fitting-view.ts#L1621-L1622)).
- **Both run branches** ([fitting-view.ts:1626-1675](function-views/src/fitting-view.ts#L1626-L1675)) pass the selected `method` and its settings map:
  `getFittedParamsFinalized({... method, settings ...})` (Diff Studio) and `runOptimizerFinalized({... method, settings ...})` (standard path + fallback).
- Generalize the `Issues` button switch ([fitting-view.ts:644-653](function-views/src/fitting-view.ts#L644-L653)).
- Remove the Gradient Descent dead code (`gradDescentSettings`, `generateGradDescentSettingsInputs`, its validator).

### 5.9 Diff Studio (IVP models) — in scope
Diff Studio models go through a separate branch ([fitting-view.ts:1626-1662](function-views/src/fitting-view.ts#L1626-L1662)):
`getFittedParamsFinalized` → `getFittedParams` (spins up a `./workers/basic` worker pool, splits starting points into batches)
→ worker `basic.ts` → `fit(task, start)`, which builds a `costFunc` (solves the IVP via `applyPipeline`, metric `mad`/`rmse`)
and calls the **async** `optimizeNM` ([diff-studio/fitting-utils.ts:190](function-views/src/fitting/diff-studio/fitting-utils.ts#L190)).
On a worker error it falls back to `runOptimizerFinalized` (main-arm), which already supports the methods.

Because this uses the async optimizer inside the worker, **codegen / sync twins are not needed** — the sci-comp
adapter (`minimizeAsync`) plugs in directly. Changes:
- `diff-studio/defs.ts`: add `method: string` to `NelderMeadInput` (a serializable `METHOD` value).
- `diff-studio/nelder-mead.ts`: add a `method` parameter to `getFittedParams`/`getFittedParamsFinalized`, put it into `task`.
- `diff-studio/fitting-utils.ts` (`fit`): instead of the hard `optimizeNM` — take the optimizer from the registry by `task.method`;
  when `wantsBounds`, build `OptimizerBounds` from `task.bounds`+`task.variedInputNames` (in the same order as vector `x`)
  and pass it to the optimizer. The internal `costOutside` cutoff remains as a safety net over the native/penalty bounds.
- `diff-studio/workers/basic.ts`: logic unchanged — `fit()` reads `task.method` itself.
- **Requirement**: `optimizer-registry.ts` (and everything it pulls in: the adapter, sci-comp, the descriptors) must stay
  **worker-safe** — no `grok.*`/`ui.*`/runtime `DG` (type-only DG imports only). sci-comp and the adapter are pure TS, so this is achievable.

### 5.10 Out of scope (explicit)
- The **standard path's worker pool** (`fitting/worker/`, `fitting.worker.ts`) uses the **sync** `optimizeNMSync` (codegen).
  Adding new methods there would require sync twins (`npm run update-codegen`) + settings serialization — a separate follow-up.
  In this phase the new methods on the standard path run main-arm only; the NM codegen is untouched.

## 6. Files affected

| File | Change |
|---|---|
| `package.json` | + dependency `@datagrok-libraries/sci-comp` |
| `.../fitting/optimizer-misc.ts` | + `OptimizerBounds`, extend `IOptimizer` (optional `bounds`) |
| `.../fitting/optimizer-sci-comp.ts` | **new**: adapter + descriptors + per-method `buildSettings` |
| `.../fitting/optimizer-registry.ts` | **new**: `OPTIMIZERS: Map<METHOD, OptimizerDescriptor>` |
| `.../fitting/constants.ts` | METHOD +4, −`GRAD_DESC`, `methodTooltip` |
| `.../fitting/optimizer-sampler.ts` (or nearby) | + helper for `OptimizerBounds` from `getAccData` |
| `.../fitting/optimizer.ts` | generalize the dispatcher `performOptimization` |
| `.../fitting/worker/executor.ts` | `ExecutorArgs`+optimizer/bounds; `MainExecutor.run`; `canHandle` gate |
| `.../fitting/optimizer-api.ts` | `OptimizerParams.method`; `resolveParams`/`runRawOptimizer` via the registry |
| `.../fitting/diff-studio/defs.ts` | `NelderMeadInput.method: string` |
| `.../fitting/diff-studio/nelder-mead.ts` | `getFittedParams(Finalized)` += `method`, put into task |
| `.../fitting/diff-studio/fitting-utils.ts` | `fit()` — optimizer from the registry by `task.method` + bounds |
| `.../fitting/diff-studio/workers/basic.ts` | no logic change (fit reads task.method) |
| `.../fitting-view.ts` | dropdown, settings generation/validation, dispatch (both branches + `method`), GRAD_DESC cleanup |
| `.../compute-utils/CHANGELOG.md` (if maintained) | `v.next` entry |

## 7. Testing and verification

1. **Build**: in `libraries/compute-utils` — `grok build` (Turborepo builds sci-comp first), then `pnpm run lint`.
2. **Adapter unit tests** (next to the existing fitting tests in compute-utils): for each new method,
   fit the parameters of a known analytic model (e.g. exponential decay or a quadratic with a known minimum)
   and check: `cost` below a threshold, `point` near the true value, `iterCosts`/`iterCount` populated;
   check the `undefined→penalty` mapping, `threshold→onIteration`, and bounds passing (the point does not
   leave the box for L-BFGS-B). The ground truth is sci-comp's own tests.
3. **Bounds check**: a test with formula bounds (`formulaBounds`) — verify that ±∞ is set for those
   dimensions and the protective undefined cutoff kicks in.
4. **Diff Studio**: an integration test on a known IVP model — fit the parameters with each method via
   `getFittedParamsFinalized` (in the worker); check convergence and that the registry loads correctly in
   the worker (worker-safe). Check the fallback to main-arm on a worker error.
5. **Manual UI check**: open the Fitting View (a) on a regular model and (b) on a Diff Studio model, select
   `Nelder-Mead`, `L-BFGS-B`, `PSO`, `L-BFGS`, `Adam` in turn; verify that the method's settings inputs
   change, the run completes, the "Fitting profile" chart is drawn, and the result converges.
6. **NM regression**: verify that Nelder-Mead (standard path main/worker with `//meta.workerSafe`, and the
   Diff Studio worker) is unchanged in behavior.

## 8. Resolved questions and risks

### Resolved (see the decisions table in §2)
1. **Iteration defaults** — conservative values in the descriptors (reference NM=50), NOT sci-comp's 1000/10000.
2. **Penalty μ** — a separate method setting (L-BFGS/Adam/PSO), default 1000 → `penaltyOptions.mu`.
3. **Gradient cost** — capped with a default `maxFunctionEvaluations`/`maxIterations`.
4. **Nelder-Mead** — keep the self-contained compute-utils NM; new methods via sci-comp.

### Resolved by default (documented, no special input required)
5. **Settings subset in the UI**: start with a minimally useful set of fields per method, expand as needed.
6. **`converged` flag**: sci-comp returns it, but `Extremum` doesn't store it — ignored for now (the `Extremum` shape is unchanged).
7. **Cancellation**: forward the cancel check into the adapter's `onIteration` (the way the standard path checks `pi.canceled` between seeds).
8. **Units**: `costHistory` (`Float64Array`) → `iterCosts` (`number[]`) — converted in the adapter.

### Implementation risks (verified in code)
9. **Registry worker-safety**: `optimizer-registry.ts` (+ adapter, sci-comp, jstat) must bundle and run
   correctly inside the Diff Studio worker without main-thread globals. sci-comp is pure TS; check the worker
   bundle size (it pulls in `jstat`).
10. **Bounds order consistency** in Diff Studio: build `OptimizerBounds` in the same order as vector `x`
    (`variedInputNames`/`boundsIdx`), as in `sampleParamsWithFormulaBounds`.

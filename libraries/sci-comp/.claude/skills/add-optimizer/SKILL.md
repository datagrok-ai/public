---
name: add-optimizer
description: Add a new single-objective optimization algorithm to the sci-comp library. Use when the user asks to implement a new optimizer (e.g. gradient descent, BFGS, simulated annealing, differential evolution).
user-invocable: true
---

# Add a new single-objective optimizer

The user wants to add a new optimizer: **$ARGUMENTS**

Follow these steps exactly. Do not skip any step.

## Step 1: Understand the algorithm

Before writing code, research the algorithm:
- What are its tunable hyperparameters?
- Is it population-based or single-point?
- Does it require gradients?
- What is the iteration logic?

## Step 2: Read the reference implementation

Read these files to understand codebase patterns:
- `src/optimization/single-objective/optimizers/nelder-mead.ts` — reference optimizer
- `src/optimization/single-objective/__tests__/nelder-mead.test.ts` — reference tests (the new optimizer MUST include all the same test cases)
- `src/optimization/single-objective/__tests__/helpers.ts` — test functions and helpers
- `src/optimization/single-objective/examples/unconstrained.ts` — reference example

## Step 3: Create the optimizer file

Create `src/optimization/single-objective/optimizers/<name>.ts` following the pattern from `nelder-mead.ts`.

### Rules (beyond what's visible in the reference)

- `runInternal` receives the already-penalized objective — do NOT handle constraints
- `runInternal` receives settings with defaults already applied via `withDefaults`
- Set `converged = true` only if a convergence criterion was met, not just maxIterations
- **Pre-allocate all buffers before the loop** — zero allocations inside the iteration body
- Helper methods must write into caller-provided buffers (out-parameter pattern), not allocate new arrays
- `costHistory`: pre-allocate `Float64Array(maxIter)` with a separate `costLen` counter; return `costHistory.subarray(0, costLen)`

## Step 4: Export from index.ts

Edit `src/optimization/single-objective/index.ts`:

1. Add export for the class and settings type (in the "Built-in optimizers" section):
```typescript
export {<Name>} from './optimizers/<name>';
export type {<Name>Settings} from './optimizers/<name>';
```

2. Add auto-registration (in the "Auto-register" section at the bottom):
```typescript
import {<Name>} from './optimizers/<name>';
registerOptimizer('<kebab-name>', () => new <Name>());
```

## Step 5: Write tests

Create `src/optimization/single-objective/__tests__/<name>.test.ts`.

The new test file MUST replicate ALL test cases from `nelder-mead.test.ts` — every `describe` block, each with both `sync` and `async` variants.

You may adjust:
- `x0` starting points (if the algorithm needs a closer start)
- `maxIterations`, `tolerance`, and algorithm-specific settings
- Precision in `toBeCloseTo` / `expectPointClose` (if the algorithm is less precise)

You must NOT:
- Remove any test group
- Change expected values or expected points

## Step 6: Create example file

Create `src/optimization/single-objective/examples/<name>.ts` following the structure of `unconstrained.ts`. Must show `minimize`, `maximize`, and at least one constrained example with `boxConstraints` + `applyPenalty`.

## Step 7: Verify

Run in order:
1. `npm run lint-fix`
2. `npm run build` — must compile without errors
3. `npm test` — run all tests

**CRITICAL:** If any tests fail, do NOT silently fix or skip them. Instead:
1. Collect the full list of failing test names and reasons
2. Present the list to the user
3. **Wait for the user's response** — do NOT proceed until the user explicitly approves a course of action
This rule applies to every test run, including re-runs after fixes.

## Step 8: Add to benchmarks

Add the new optimizer to `src/optimization/single-objective/benchmarks/unconstrained-benchmarks.ts`:

1. Import the new optimizer class
2. Add an entry to the `optimizers` array:
```typescript
{
  name: '<Name>',
  optimizer: new <Name>(),
  settings: {maxIterations: 10_000, /* algorithm-specific defaults */},
},
```
3. Run: `npx tsx src/optimization/single-objective/benchmarks/unconstrained-benchmarks.ts`
4. Save results to `src/optimization/single-objective/benchmarks/unconstrained-benchmarks.md` (reformat as markdown tables, see existing file)
5. Update the problem count in the banner if it changed

## Step 9: Update CLAUDE.md

In `CLAUDE.md`, update the architecture tree — add the new optimizer file under `optimizers/`.

## Step 10: Add to README.md

In `README.md`, add the new optimizer to the list under "Single-objective" section, with a Wikipedia or reference link.

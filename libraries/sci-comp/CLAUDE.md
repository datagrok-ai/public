# CLAUDE.md

`@datagrok-libraries/sci-comp` is a pure TypeScript library of numerical methods for the Datagrok platform. It has **no dependency on datagrok-api**.

## Domains

```
index.ts                      # Re-exports {singleObjective, multiObjective, timeSeries, stats, nca} namespaces
src/optimization/             # Single- and multi-objective solvers — see src/optimization/CLAUDE.md
src/stats/                    # Statistical tests, distributions, multiple-comparison — see src/stats/CLAUDE.md
src/time-series/              # Time-series feature extraction — see src/time-series/CLAUDE.md
src/nca/                      # Non-compartmental analysis (PK) — see src/nca/CLAUDE.md
```

Domain-specific architecture trees, design patterns, and conventions live in the per-domain `CLAUDE.md`. Open the relevant one when working in that subtree — Claude Code loads nested `CLAUDE.md` files automatically.

## Cross-cutting invariants

- **No `datagrok-api` dependency**: this library stays platform-agnostic. DataFrame adapters, UI, worker pools belong in the consuming package (e.g. `packages/NCA/`), not here.
- **`Float64Array` for vectors**: all numeric point vectors use `Float64Array`, never `number[]`. Internal accumulators also stay in `Float64Array`.
- **Namespace re-exports** (`index.ts`): the public API exposes one namespace per domain to avoid name collisions. Consumers import as `import {singleObjective, stats, timeSeries, nca} from '@datagrok-libraries/sci-comp'`. New domains follow the same pattern.

## Adding a new method

Follow [`.claude/rules/new-method-checklist.md`](.claude/rules/new-method-checklist.md) — every item is required:

1. Tests in the relevant `__tests__/` (sync + async where applicable).
2. TSDoc on every public type, interface, function, class.
3. Section README in the same directory (or nearest parent).
4. Top-level `README.md`.
5. **The CLAUDE.md for the relevant domain** — update the architecture tree.
6. `npm run lint-fix && npm run build && npm test` all pass.

## Testing

`npm test` runs jest (ts-jest transform) over `src/**/__tests__/`, per [`jest.config.js`](jest.config.js).

```bash
npm test                      # whole suite (48 suites, ~1270 tests, ~20 s)
npx jest src/nca              # one domain — path substring
npx jest -t 'lambda_z'        # one test or describe block — name substring
npx jest -u                   # update snapshots: see src/nca/CLAUDE.md before doing this
npx jest --coverage
```

This library pins its own `typescript` (`npm:typescript@^5.9.3`) instead of taking `catalog:`, which the
workspace catalog holds at `^7.0.2`. TypeScript 7's CJS entry point is a version stub — the compiler API
moved to `./unstable/*` — and ts-jest 29 requires the TS 5 API (peer range `>=4.3 <6`). Without the pin every
suite dies in the transform with `Cannot read properties of undefined (reading 'Node16')`.

The pin is confined to the test path: `grok tsc` resolves TypeScript from `@datagrok/build-config`'s own
directory, not from this package, so `npm run build` still compiles with the catalog's TypeScript 7.

## Skills

- `/add-optimizer <algorithm name>` — scaffold a new single-objective optimizer (see [`src/optimization/CLAUDE.md`](src/optimization/CLAUDE.md)).

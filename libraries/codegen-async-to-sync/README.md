# @datagrok-libraries/codegen-async-to-sync

Generates a sync mirror of an async TypeScript source file via a `ts-morph` AST
transform. The async file is the source of truth; the sync sibling is regenerated
on demand and checked into git alongside it.

This is a **dev-only** tool. Add it to a library's `devDependencies` and invoke
the bundled CLI from `package.json` scripts.

## Install

```bash
npm install --save-dev @datagrok-libraries/codegen-async-to-sync
```

```jsonc
// package.json
"scripts": {
  "update-codegen": "codegen-async-to-sync --roots src",
  "codegen:check":  "codegen-async-to-sync --check --roots src"
}
```

`--roots` accepts one or more directories to scan (relative to CWD). Files are
only considered if their leading comment block contains an `@async-source`
directive.

`--check` mode regenerates in memory and compares with the committed sync sibling.
Exits non-zero if any sync file would change — for CI drift detection.

## Directive grammar

In the **leading comment block** of an async source file:

```ts
// @async-source: <output-filename>      // required — names the generated sibling
// @codegen-rename: <old>=<new>          // optional, repeatable
```

`@codegen-rename` renames top-level function and variable declarations AND
**named-import bindings**. The latter is how you route an async import
(`runLineSearchAsync`) to its sync sibling (`runLineSearch`) in the generated
output — the rename is applied to the import specifier and the rebuild step
emits the renamed name, assuming the sibling module exports it (standard pattern
when both sides follow the async/sync-twin convention).

The codegen rejects two malformed renames at parse time:
- `foo=foo` — source equals target (no-op typo).
- Two renames with the same target (`a=shared`, `b=shared`) — ambiguous.

Inline directives in the body:

```ts
const debugOnly = true;            // @async-only         — drops just this line

// @async-only-begin
if (asyncOnlyFlag) await yield();  // these lines (and the markers)
await flushQueue();                // are all dropped
// @async-only-end
```

Block form is for async-only logic that spans multiple statements — the
single-line form would strip only the head and leave a syntactically broken
tail. Unmatched `-begin` strips to EOF (intentional — surfaces as a parse
error downstream); unmatched `-end` is a no-op.

## Transformation rules

- `async function`, `async function expression`, and `async () => …` arrows
  — **top-level OR nested** — have their `async` modifier stripped.
- `await expr` → `expr` (anywhere in the source, including inside nested closures).
- `Promise<T>` → `T` (iterates to a fixed point for nested cases).
- `@codegen-rename` applies to declared identifiers, their references, and
  named-import bindings (above).
- Top-level declarations that are **not async** — functions, variables,
  **interfaces, type aliases** — are not copied to the sync output; if the
  copied code references them they are imported back from the source file via
  `import {sibling} from './<srcStem>'`.
- `import * as X from 'mod'` — passed through unchanged when the binding is
  referenced in the sync output, dropped when unused. Namespace bindings are
  not subject to `@codegen-rename`.
- Generated file starts with `/* eslint-disable */` — generated output
  shouldn't be lint-gated.

## Consumer patterns

### File-to-file (top-level async functions)

The async file declares top-level `async function` or `async () => ...` decls.
The generated sibling exports the sync equivalents.

```ts
// foo.ts (source of truth)
// @async-source: foo-sync.ts
export async function fooAsync(x: number): Promise<number> { return await bar(x); }

// foo-sync.ts (GENERATED — do not edit)
```

### Class shim → top-level driver

For class methods, refactor the algorithm body into a top-level function,
have the class delegate to it, and codegen the sync top-level mirror.

```ts
// optimizer.ts
export class MyOptimizer {
  protected runInternal(...)            { return runMyOptimizerSync(...);  }
  protected async runInternalAsync(...) { return runMyOptimizerAsync(...); }
}

// optimizer-driver.ts (source of truth)
// @async-source: optimizer-driver-sync.ts
// @codegen-rename: runMyOptimizerAsync=runMyOptimizerSync
export async function runMyOptimizerAsync(...) { ... }

// optimizer-driver-sync.ts (GENERATED)
```

This is the LBFGS-B pattern in `@datagrok-libraries/sci-comp`.

## Known limitations

The codegen does **not** currently handle the following — outputs may be
silently broken (S) or produce a clear error / clean failure (E):

| Pattern | Status | Workaround |
|---|---|---|
| `export default async function …` | E — throws "no async top-level decls" | Use named exports |
| `async function*` (async generators) | **S** — `async` stripped but generator stays; yields `Promise<T>` instead of `T` | Refactor; await values explicitly inside a regular async fn |
| `for await (const x of …)` | E — emits invalid syntax (`for await` in non-async fn) | Use a regular `for-of` if the iterable is sync, otherwise refactor |
| `const a = sync, b = async …` (multi-declarator mixed async) | **S** — treated as fully async, sibling-discovery breaks | Split into separate statements |
| `import {x as y} from 'm'` + `@codegen-rename: y=z` | E — import dropped from output (rebuild filter doesn't track local bindings) | Don't combine alias with rename; pick one |
| `export {x} from './y'` (re-exports) | **S** — dropped from sync output | Convert to `import` + `export` |
| Sync factory returning async closure (`function f() { return async () => …; }`) | not supported as `@async-source` — the codegen requires top-level **async** | Hand-maintain the sync twin |

These are tracked for follow-up; open an issue if you hit one.

## Exit codes

- `0` — success (or no codegen sources found).
- `1` — `--check` mode detected drift; sync mirror differs from what would be generated.
- `2` — invocation error (bad flag, source file invalid, etc).

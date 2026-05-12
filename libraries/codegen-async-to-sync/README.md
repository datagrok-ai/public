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
// @codegen-rename: <old>=<new>          // optional, repeatable — renames identifiers
```

Inline directive (any line):

```ts
const debugOnly = true; // @async-only   // line is stripped from sync output
```

## Transformation rules

- `async function` / `async () => ...` → `async` removed.
- `await expr` → `expr`.
- `Promise<T>` → `T` (iterates to a fixed point for nested cases).
- Apply `@codegen-rename` to declared identifiers and their references.
- Top-level declarations that are **not async** are not copied to the sync output;
  if the copied code references them they are imported from the source file
  via `import {sibling} from './<srcStem>'`.

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
  protected runInternal(...)      { return runMyOptimizerSync(...);  }
  protected async runInternalAsync(...) { return runMyOptimizerAsync(...); }
}

// optimizer-driver.ts (source of truth)
// @async-source: optimizer-driver-sync.ts
export async function runMyOptimizerAsync(...) { ... }

// optimizer-driver-sync.ts (GENERATED)
```

This is the LBFGS-B pattern in `@datagrok-libraries/sci-comp`.

## Exit codes

- `0` — success (or no codegen sources found).
- `1` — `--check` mode detected drift; sync mirror differs from what would be generated.
- `2` — invocation error (bad flag, source file invalid, etc).

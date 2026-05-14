---
name: datagrok-calc-column
description: Add a calculated, formula-based column to a dataframe inside a datagrok-exec block. Use whenever the user asks to compute, derive, add, or create a new column from existing columns — LipE, ratios, log/round, heavy atom count, any expression in the Datagrok formula DSL. Replaces hand-written addNewFloat/addNewInt + for-loop with a single formula-attached column that recomputes when source columns change.
---

# datagrok-calc-column

Use `grokky.addCalculatedColumn(...)` inside a `datagrok-exec` block to add a
formula-driven column. The formula stays attached to the column, so edits to
source columns trigger automatic recompute.

Do **not** use `t.columns.addNewFloat('LipE')` plus a manual `for` loop — that
produces a static column with no formula and no recompute.

## Signature

```ts
await grokky.addCalculatedColumn(
  df,                // DG.DataFrame  (use `t` in a TableView)
  name,              // string
  formula,           // string — DSL, see below
  type?              // 'auto' (default) | 'double' | 'int' | 'string' | 'bool' | 'datetime'
): Promise<DG.Column>
```

## Formula DSL — quick reference

| Form                 | Meaning                                            |
|----------------------|----------------------------------------------------|
| `${ColumnName}`      | per-row scalar reference, exact case               |
| `$[ColumnName]`      | whole-column reference (use inside aggregates)     |
| `row`                | current row index                                  |
| `+ - * / ^`          | operators                                          |
| `Round`, `RoundFloat`, `Log`, `Log10`, `Sqrt`, `Abs` | math functions      |
| `Avg($[col])`, `Sum($[col])`, `Min`, `Max`, `StDev` | aggregates over a whole column |
| `Chem:getInchis(${molecule})`, `HeavyAtomCount(${molecule})` | only functions tagged `vectorFunc: 'true'` work here — see [Which functions work in a formula](#which-functions-work-in-a-formula) |
| `"text"`             | string literal (double quotes)                     |

Both operator and function forms parse: `${pIC50} - ${cLogP}` ≡ `Sub(${pIC50}, ${cLogP})`.

Full function catalog: `help/transform/add-new-column.md`.

## Examples

```datagrok-exec
// Heavy atom count (requires Chem package)
await grokky.addCalculatedColumn(t, 'HAC', 'HeavyAtomCount(${molecule})', 'int');
```

```datagrok-exec
// Log of activity
await grokky.addCalculatedColumn(t, 'logIC50', 'Log10(${IC50})');
```

## Which functions work in a formula

A package function (`Chem:foo`, `Admetica:bar`, ...) can be used inside an
`addCalculatedColumn` formula only if it is tagged `vectorFunc: 'true'` in its
metadata. Everything else must be called imperatively via
`grok.functions.call(...)` in a separate `datagrok-exec` block.

Check the tag before emitting the formula — for example, with the MCP
`get_function` call, or by reading the package's `package.ts`. Do not infer it
from the function name.

## What this skill does NOT do

- It does not filter, sort, or select rows — those are separate skills.
- It does not add viewers — use `grok.shell.addViewer(...)` or a future viewer skill.
- It does not run server-side scripts — use `grok.functions.call(...)` or the script's registered name.

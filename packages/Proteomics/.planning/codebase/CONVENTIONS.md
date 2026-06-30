# Coding Conventions

**Analysis Date:** 2026-05-11

## Code Style

- **Indentation:** 2 spaces
- **Quotes:** Single quotes (`'`) for strings; backticks for template literals
- **Semicolons:** Required at end of every statement
- **Line endings:** CRLF (Windows-style) — repository-wide convention
- **Line length:** 120-char soft limit (ESLint, Google style)
- **TypeScript:** Strict mode, ES2020 target, decorators enabled

Enforced by `.eslintrc.json` (Google base config) and `tsconfig.json`.

## File Naming

- **Source files:** kebab-case — `differential-expression.ts`, `maxquant-parser.ts`, `qc-dashboard.ts`
- **No file-type suffix:** not `maxquant-parser.service.ts`, just `maxquant-parser.ts`
- **R scripts:** snake_case ending in `.R` — `limma_de.R`, `vsn_normalize.R`
- **Test files:** match the file they test — `src/tests/qc-dashboard.ts` for `src/viewers/qc-dashboard.ts`

## Function Naming

The package uses a deliberate three-way split for the most common function shapes:

| Prefix | Returns | Side effect | Example |
|--------|---------|-------------|---------|
| `showX(df)` | `void` (opens a `ui.dialog(...)`) | Opens a modal dialog; mutates `df` on OK | `showDEDialog`, `showNormalizationDialog`, `showAnnotationDialog` |
| `openX(df)` | `void` (creates a view) | Adds a table view or dashboard to the shell | `openQcDashboard`, `openEnrichmentVisualization` |
| `createX(df, opts?)` | `DG.Viewer` | Pure factory — no shell side effects | `createVolcanoPlot`, `createExpressionHeatmap`, `createPcaPlot` |

Other patterns:
- **Parsers:** `parseXText(text)` (synchronous, takes raw string) and `parseX(file)` (async wrapper)
- **Computations:** verb + noun — `medianNormalize`, `imputeKnn`, `runDifferentialExpression`, `computePCA`
- **Helpers:** `findColumn`, `getGroups`, `setGroups` — small camelCase verbs

## Class and Constant Naming

- **Classes:** PascalCase — `PackageFunctions`, `ProteomicsPackageDetectors`
- **Interfaces:** PascalCase — `GroupAssignment`, `EnrichmentOptions`
- **Constants:** SCREAMING_SNAKE_CASE inside an `as const` object — single-source-of-truth pattern:
  ```ts
  export const SEMTYPE = {
    PROTEIN_ID: 'Proteomics-ProteinId',
    GENE_SYMBOL: 'Proteomics-GeneSymbol',
    LOG2FC: 'Proteomics-Log2FC',
    P_VALUE: 'Proteomics-PValue',
    INTENSITY: 'Proteomics-Intensity',
  } as const;
  ```

## Function Registration — Metadata Comments

Every platform-exposed function lives as a static method on `PackageFunctions` in `src/package.ts`, decorated with `@grok.decorators.func({...})`. Top-menu entries follow a strict path syntax:

```ts
@grok.decorators.func({'top-menu': 'Proteomics | Analyze | Normalize...'})
static async normalizeProteomics(): Promise<void> {
  const df = grok.shell.tv?.dataFrame;
  if (!df) { grok.shell.warning('No table open'); return; }
  showNormalizationDialog(df);
}
```

The trailing `...` in the menu path is convention for items that open a dialog (vs. immediate-action items). The `grok api` build step parses these decorators to generate `src/package.g.ts` and `src/package-api.ts`.

## DataFrame Tags — Workflow State

State accumulates as tags on the source `DG.DataFrame`. All package tags use the `proteomics.` prefix:

| Tag | Value | Set by |
|-----|-------|--------|
| `proteomics.source` | `'maxquant' \| 'spectronaut' \| 'generic'` | Parsers on import |
| `proteomics.groups` | JSON-encoded `GroupAssignment` | `setGroups(df, groups)` after Annotate Experiment |
| `proteomics.normalized` | `'true'` | Any normalization completion |
| `proteomics.preNormalized` | `'true'` | Spectronaut import (if input was already normalized) |
| `proteomics.imputed` | `'true'` | Any imputation completion |
| `proteomics.de_complete` | `'true'` | DE pipeline completion |
| `proteomics.de_method` | `'limma' \| 'deqms' \| 't-test'` | DE on completion (which method actually ran) |
| `proteomics.enrichment` | `'true'` | Enrichment DataFrame marker |

Downstream functions check these tags as preconditions (e.g., heatmap reads `de_complete` before rendering).

## User Feedback Patterns

The shell-level notification functions are used by category:

- **`grok.shell.warning(msg)`** (~35 uses) — Precondition violations, missing inputs, fallback notifications ("R unavailable, using client-side t-test")
- **`grok.shell.info(msg)`** (~14 uses) — Successful completion ("Normalization complete, 1,247 proteins processed")
- **`grok.shell.error(msg)`** (~3 uses) — Hard failures the user must act on
- **`console.warn(msg)`** — Recoverable internal errors the user doesn't need to see directly (e.g., Dendrogram fallback)

## Input Tooltips

Every dialog input gets a `.setTooltip(...)` call immediately after creation. Pattern from `src/analysis/differential-expression.ts`:

```ts
const methodInput = ui.input.choice('Method', {value: 'limma', items: ['limma', 'DEqMS', 't-test']});
methodInput.setTooltip('limma: moderated t-test; DEqMS: peptide-count-weighted variance; t-test: client-side Welch\'s t-test');
```

(Higher-level dialog/viewer help text is largely missing — see CONCERNS.md "Analytics-Level Help Missing".)

## Progress Indicators

All operations expected to take >100ms wrap their work in `DG.TaskBarProgressIndicator`:

```ts
const pi = DG.TaskBarProgressIndicator.create('Computing differential expression...');
try {
  pi.update(20, 'Building expression matrix...');
  const exprDf = buildExpressionDf(df, groups);
  pi.update(60, 'Running limma...');
  await grok.functions.call('Proteomics:LimmaDE', {...});
  pi.update(100, 'Done');
} finally {
  pi.close();
}
```

The `try/finally` is required so the indicator closes even if the operation throws. Mid-step `pi.update(percent, message)` calls are used inconsistently — see CONCERNS.md.

## Column Lookup Strategy

The package never grabs columns by name directly. Instead, it uses the two-tier strategy in `src/utils/column-detection.ts`:

1. **Semantic type first:** `df.columns.bySemType(SEMTYPE.GENE_SYMBOL)`
2. **Name-hint fallback:** match against a list of candidate names if semantic type isn't tagged

```ts
const geneCol = findColumn(df, SEMTYPE.GENE_SYMBOL, ['gene name', 'gene symbol', 'genes']);
if (!geneCol) { grok.shell.warning('No gene symbol column found'); return; }
```

This makes the package robust to imports from new vendors that haven't been semantically typed yet.

## Error Handling

Cascading try/catch with graceful degradation is the dominant pattern:

```ts
try {
  await grok.functions.call('Proteomics:LimmaDE', {...});  // server-side R
} catch (e) {
  grok.shell.warning('R limma unavailable, falling back to client-side t-test');
  return runDifferentialExpression(df, groups);  // pure JS
}
```

Almost every R-dependent path has a JS fallback. The DE pipeline has a three-level chain: DEqMS → limma → client-side t-test.

## File Organization

Module structure (each in its own file under `src/`):

- **`src/analysis/`** — Pure computation that mutates the DataFrame; one file per analytical step
- **`src/parsers/`** — File-format → DataFrame; one file per vendor + a `shared-utils.ts` for common helpers
- **`src/viewers/`** — `DG.Viewer` factories; one file per viewer or dashboard
- **`src/panels/`** — Context-panel widgets
- **`src/utils/`** — Cross-cutting constants and helpers (semantic types, column lookup)
- **`src/tests/`** — Test suites; one file per domain area

No barrel files (no `index.ts` re-exports). Imports go directly to the source file.

## Imports

- **Platform namespaces** — three required imports at the top of every TS source:
  ```ts
  import * as grok from 'datagrok-api/grok';
  import * as ui from 'datagrok-api/ui';
  import * as DG from 'datagrok-api/dg';
  ```
- **Internal modules** — relative paths within the package (no path aliases configured)
- **No default exports** — all exports are named

## Comments

Sparse. The codebase generally favors well-named identifiers over inline comments. The conventions for the few comments that exist:

- **JSDoc-style** for public exported functions when their behavior is non-obvious
- **`// TODO:` / `// FIXME:`** are rare; when they exist they should be promoted to `.planning/todos/`
- **`//name:` / `//top-menu:` / `//description:`** metadata comments at top of R scripts and (historically) above some functions — now superseded by `@grok.decorators.*` for TS

## Things to Avoid

- **Editing `src/package.g.ts` or `src/package-api.ts`** — auto-generated; overwritten by `grok api`
- **Importing `datagrok-api/dg`, `grok`, `ui`, `rxjs` as bundled dependencies** — they are webpack externals provided by the platform
- **Hand-coding column names** — use `findColumn(df, SEMTYPE.X, [...])` instead
- **Bare `fetch()` for external URLs** — should use `grok.dapi.fetchProxy()` (currently violated in UniProt and g:Profiler paths — see CONCERNS.md)
- **Mutating a DataFrame in a way that breaks downstream preconditions** — set the corresponding `proteomics.*` tag so downstream checks can react

---

*Conventions analysis: 2026-05-11*

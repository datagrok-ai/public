# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

Datagrok extension for mass-spectrometry proteomics analysis. Imports search-engine output
(MaxQuant, Spectronaut, FragPipe, generic CSV/TSV), runs the analysis pipeline
(annotate groups -> normalize -> impute -> differential expression -> enrichment),
and adds the matching viewers (volcano, heatmap, PCA, QC dashboard, enrichment charts).
Exposed in the platform under the top-menu **Proteomics**.

Auth / credentials: none. The package only calls public REST APIs (UniProt, g:Profiler).
R scripts (`scripts/limma_de.R`, `scripts/deqms_de.R`, `scripts/vsn_normalize.R`) run server-side
via `grok.functions.call(...)` and need a configured Datagrok R compute environment; every
R path has a client-side TypeScript fallback.

## Architecture

Functional pipeline. Each step mutates a shared `DG.DataFrame` in-place and records its
completion as a `proteomics.*` tag on the DataFrame. Downstream functions read those tags
as preconditions. There is no separate state store.

- **`src/package.ts`** — entry point. `PackageFunctions` class with `@grok.decorators.func`
  / `init` / `panel` methods. Each handler validates that a table is open, checks any
  required `proteomics.*` precondition tag, then delegates to `analysis/` or `viewers/`.
  The **Proteomics** menu is built two ways (see "Menu architecture" below).
- **`src/menu.ts`** — builds the full **Proteomics** ribbon menu (Import + Annotate /
  Analyze / Visualize / Share) on a table view's `ribbonMenu`, with the `isEnabled` grey-out:
  sample-level items disable themselves on a Spectronaut Candidates table (no per-sample
  intensities). The `buildProteomicsTopMenu` autostart in `package.ts` calls
  `buildProteomicsRibbonMenu` once per Proteomics table view.

### Menu architecture (read before touching the menu — non-obvious, verified live)

The **Proteomics** menu lives in two independent hosts. A decorator `top-menu` populates
`grok.shell.topMenu` (the **left-bar main menu**) ONLY; it does NOT populate a table view's
ribbon. So:

- **Left-bar main menu** = decorator `top-menu` strings on the `importX` handlers in
  `package.ts`, e.g. `Proteomics | Import | Spectronaut Candidates...`. This is the
  always-available entry point — reachable with **no table open** (the ribbon shows no
  package menu on the Home/start view). Declaration order = menu order (Candidates first,
  then Report). Only Import is here, because import is the only thing you do without a table.
- **Ribbon (table view)** = built programmatically in `menu.ts` onto `view.ribbonMenu`, and
  only when a Proteomics table (`proteomics.source` tag) is active. It builds the WHOLE menu
  including its own Import group — the decorator Import does not reach the ribbon, and a
  view's own `ribbonMenu` 'Proteomics' group is what the ribbon shows. Built once per view;
  do NOT rebuild on view activation (`isEnabled` re-evaluates on open — rebuilding the menu
  duplicates it, which was a real bug).

To add a menu item: if it's an importer, add a decorator `top-menu` handler in `package.ts`
(left-bar menu) AND an `imp.item(...)` line in `menu.ts` (ribbon). Otherwise add the handler
in `package.ts` and an item in the right group in `menu.ts` (pass `sampleOnly` as the
`isEnabled` option if it needs per-sample data). Don't give the analysis groups decorator
`top-menu` strings — they'd show in the left-bar menu with no table context and can't grey
out.
- **`src/analysis/`** — pure computation that mutates the DataFrame and sets a completion
  tag. One file per step: `experiment-setup.ts` (groups + `getGroups`/`setGroups`),
  `normalization.ts`, `imputation.ts`, `differential-expression.ts`, `pca.ts`, `enrichment.ts`.
- **`src/parsers/`** — vendor file -> filtered, log2-transformed, semantic-typed DataFrame.
  `maxquant-parser.ts`, `spectronaut-parser.ts` (PG report — per-sample quant),
  `spectronaut-candidates-parser.ts` (Candidates report — pre-computed DE; emits the
  canonical `log2FC` / `p-value` / `adj.p-value` / `significant` shape and sets
  `proteomics.de_complete`, so the pipeline skips Annotate/Normalize/Impute/DE),
  `fragpipe-parser.ts`, `generic-parser.ts`, plus `shared-utils.ts` (log2 transform,
  log2-status detection, delimiter detection, primary-id column extraction). Parsers return
  unnamed DataFrames; the import handler in `package.ts` sets `df.name` from the filename.
- **`src/viewers/`** — `DG.Viewer` factories that read an already-analysed DataFrame.
  `volcano.ts`, `heatmap.ts`, `pca-plot.ts`, `qc-dashboard.ts` (+ `qc-computations.ts` for
  MA / CV / missingness math), `enrichment-viewers.ts`.
- **`src/panels/uniprot-panel.ts`** — context-panel widget that fetches per-protein metadata
  on cell click. Registered via `@grok.decorators.panel` with a `semType` filter.
- **`src/utils/proteomics-types.ts`** — `SEMTYPE` constants. **Single source of truth** for
  every `Proteomics-*` semantic type string. Use these; never inline the strings.
- **`src/utils/column-detection.ts`** — `findColumn(df, semType, nameHints)` and
  `findProteomicsColumns(df)`. **Use these instead of `df.col('Gene names')`** so the
  package keeps working on tables from new vendors that haven't been semantically tagged.
- **`detectors.js`** (package root, plain JS — not under `src/`) — semantic-type detectors
  loaded by the platform before the webpack bundle. Mirrors the `SEMTYPE` constants. If you
  add a new semantic type, add a detector here too.

### Pipeline shape (read before editing)

```
parse* (parsers/)
  -> grok.shell.addTableView(df)
  -> showAnnotationDialog (analysis/experiment-setup)  sets proteomics.groups
  -> showNormalizationDialog (analysis/normalization)  sets proteomics.normalized
  -> showImputationDialog    (analysis/imputation)     sets proteomics.imputed
  -> showDEDialog            (analysis/differential-expression) sets proteomics.de_complete (+ proteomics.de_method)
  -> createVolcanoPlot / createExpressionHeatmap / showEnrichmentDialog
```

Workflow state is **only** on the DataFrame as tags (see "Tags" below). Group assignments
are JSON-encoded in `proteomics.groups` — always go through `getGroups(df)` / `setGroups(df, g)`
in `analysis/experiment-setup.ts`; do not parse the tag yourself.

PCA is the one viewer that produces a **new** sample-level DataFrame (rows = samples, not
proteins). It must open in its own table view — never add the PCA viewer to the protein
table view. See the `showPcaPlot` handler in `package.ts` for the pattern.

### Function-naming convention (load-bearing)

The package splits the common shapes by prefix. Match the prefix when adding new code:

| Prefix | Shape | Side effect |
|---|---|---|
| `showX(df)` | opens a `ui.dialog(...)` | mutates `df` on OK |
| `openX(df)` | creates a table view / dashboard | adds to the shell |
| `createX(df, opts?)` | returns a `DG.Viewer` | none — pure factory |
| `parseXText(text)` | returns a `DG.DataFrame` | none — caller adds to shell + sets name |

Menu paths that open a dialog end with `...`; immediate-action items do not.

## Glossary

| Concept | Code | Meaning |
|---|---|---|
| Intensity column | `DG.Column` with `semType = SEMTYPE.INTENSITY` | One sample's per-protein abundance. Raw or log2-transformed; analysis steps operate on the `log2(...)` variants. |
| Group assignment | `GroupAssignment` interface (`analysis/experiment-setup.ts`); JSON in `proteomics.groups` tag | Maps two named groups (e.g. Control / Treatment) to the intensity columns belonging to each. Required input for normalization, imputation, DE, PCA, enrichment. |
| DE result | `log2FC`, `p-value`, `adj.p-value`, `significant` columns added by `runDifferentialExpression` / R scripts | Per-protein output of differential expression. Volcano/heatmap require these — they check `proteomics.de_complete`. |
| Fallback chain | `try { R } catch { client-side }` in `analysis/` | DE runs DEqMS -> limma -> client t-test; VSN -> quantile. R is optional; the package must still function without it. |
| Cross-DF link (enrichment) | `enrichDf.onCurrentRowChanged` subscription in `viewers/enrichment-viewers.ts` | Selecting a GO/KEGG term row in the enrichment DataFrame highlights the linked proteins in the protein DataFrame's volcano/heatmap. |

### DataFrame tags

All workflow state lives in tags prefixed `proteomics.`. Checking and setting these tags is
how the pipeline coordinates between steps — there is no other state store.

| Tag | Value | Set by |
|---|---|---|
| `proteomics.source` | `'maxquant'` / `'spectronaut'` / `'spectronaut-candidates'` / `'fragpipe'` / `'generic'` | parsers |
| `proteomics.groups` | JSON `GroupAssignment` | `setGroups()` after Annotate Experiment |
| `proteomics.normalized` | `'true'` | any normalization method |
| `proteomics.preNormalized` | `'true'` | Spectronaut parser (input was already normalized) |
| `proteomics.imputed` | `'true'` | any imputation method |
| `proteomics.de_complete` | `'true'` | DE pipeline completion, or Spectronaut Candidates parser (already-computed DE) |
| `proteomics.de_method` | `'limma'` / `'deqms'` / `'t-test'` / `'spectronaut'` | DE on completion (which method actually ran) |
| `proteomics.enrichment` | `'true'` | enrichment DataFrame marker |

## Conventions specific to this package

- **Look up columns via `findColumn` / `findProteomicsColumns`**, never by raw name. Imports
  from new vendors may not be semantically typed yet, so the helper falls back to name
  hints. Hard-coded `df.col('Gene names')` will silently miss valid columns.
- **Use `SEMTYPE.*` from `src/utils/proteomics-types.ts`** for every `Proteomics-*` string.
  These are the single source of truth and are mirrored in `detectors.js`.
- **Set the matching `proteomics.*` tag** when finishing an analysis step. Downstream
  viewers gate on it. If you add a new step, pick a tag name now and document it here.
- **Always provide a client-side fallback for R script calls.** The DE pipeline cascades
  DEqMS -> limma -> client-side Welch's t-test; normalization cascades VSN -> quantile.
  Never assume the platform R environment is configured.
- **Idempotency on re-run.** DE is guarded by checking `proteomics.de_complete` in the
  menu handler. Normalization/imputation are not currently guarded — if you add a step that
  adds columns, follow the `ensureFreshFloat` pattern in `viewers/qc-computations.ts`
  (remove the column if it exists, then add it) so re-running doesn't produce duplicates.
- **PCA opens in a separate table view.** Its sample-level DataFrame has a different row
  count from the protein DataFrame, so `createPcaPlot` returns `{viewer, pcaDf}` and the
  caller does `grok.shell.addTableView(pcaDf)` before docking the viewer. Do not add a PCA
  viewer to the protein table view.
- **Heatmap clones the DataFrame for filter isolation.** `createExpressionHeatmap` uses
  `df.clone(filter)` so its top-N sort/filter doesn't leak into the volcano sharing the
  same source DataFrame.
- **Adding a new semantic type:** update `SEMTYPE` in `src/utils/proteomics-types.ts` **and**
  add a detector in `detectors.js`. The README's "Semantic Type Detection" table is part of
  the public contract.
- **Adding a new parser, analysis step, or viewer:** the file layout for each is documented
  in `.planning/codebase/STRUCTURE.md` under "Where to Add New Code". Follow that contract
  (new file under the matching directory, register handler in `src/package.ts`, add
  precondition tag check, add test under `src/tests/`).
- **Deeper context:** `.planning/codebase/{ARCHITECTURE,STRUCTURE,CONVENTIONS,CONCERNS,
  INTEGRATIONS,TESTING,STACK}.md` are kept up to date and go well beyond this file. Read
  them when a task is non-trivial; do not duplicate their content here.

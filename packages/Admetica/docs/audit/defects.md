# Admetica — confirmed defects

Adversarially verified against the working tree (2026-09-09). Severity is user-facing. IDs AD-xx.

## High

- **AD-01 [resolved] — The custom tooltip hijacks every numeric column in the user's grid.**
  `src/utils/admetica-utils.ts:288-303`: `grid.onCellTooltip` gates only on "is a numeric table cell"
  and returns `true` even when the content is empty (non-Admetica column → `model` undefined →
  `getTooltipContent` returns `''`, :258-259) — returning true suppresses the platform default
  (documented in js-api `grid.md:118`; the canonical sample scopes by column, this handler doesn't).
  Installed on **every prediction run** from `applyPostProcessing` (:116) with the subscription
  discarded — N stacked handlers, permanent loss of default numeric tooltips on the whole grid.

- **AD-02 [resolved] — A per-call template silently rewires the whole session.**
  `performChemicalPropertyPredictions` overwrites the module-global `properties` from the call's
  `template` argument (`admetica-utils.ts:24, :60`); ten downstream consumers read the global
  (coloring :149, pie :253, tooltips :294, groups :320, column props :338, forms :473,
  `getModels` package.ts:47, the app…), and `setProperties` early-returns once set (:47) — the custom
  template sticks for the session, affecting other tables and the property panel. Related:
  `setProperties` canonicalizes on `items[0]` of the templates folder (:48-49) — benign with one
  shipped file, arbitrary as soon as the editor's Save writes user templates into the same folder
  (`admetica-editor.ts:78`).

- **AD-03 [resolved] — Hit Triage path crashes on unselected categories; empty selection predicts garbage.**
  `admeticaHT` spreads four `nullable: true` array params unguarded (`package.ts:66-71`); HitTriage
  replays `funcCall.inputs` verbatim with no defaulting (`libraries/statistics/src/compute-functions/
  dialog.ts:239-242`, `execution.ts:51-53`) — null → TypeError on spread. All-empty arrays are as bad:
  `models = ''` (`[]` defeats the `?? getModels()` fallback, package.ts:142) → backend
  `''.split(',') == ['']` → `find_model('')` substring-matches an arbitrary `.ckpt`
  (`dockerfiles/Admetica.py:127-133`) and returns predictions under a column named `''`.

- **AD-04 — Container: molecules with any RDKit warning are silently dropped to NaN.** *(rejected 2026-09-09: by design — any RDKit warning is intentionally treated as malformed)*
  `is_malformed` (`Admetica.py:35-58`) treats ANY captured stderr as malformed (:52) even when the
  molecule parses fine; the value becomes NaN (:86, :122-123). The `dup2` stderr-swap window is
  process-wide, so unrelated stderr during the window is misattributed too. Silent data loss that
  feeds AD-05.

## Medium

- **AD-05 [resolved] — LD50 conversion turns nulls into plausible-looking numbers.**
  `convertLD50` (`admetica-utils.ts:39`): for a null LD50 cell, `Math.pow(10, -null)` is `1`, so the
  row becomes `1 × MW × 1000` mg/kg — a realistic-looking wrong value; null MW gives NaN. Backend
  NaNs are common (AD-04). Also the only self-referential `Column.init` in the repo (reads
  `ldCol.get(i)` while re-initializing the same column) — works today, undocumented semantics.
  *(verification round: missing floats read back as the FLOAT_NULL sentinel, not `null` — the fix uses
  `Column.isNone`, and the same sentinel-as-0.000 rendering was fixed in the molecule panel.)*

- **AD-06 [resolved] — Hand-rolled CSV payload breaks on real-world values.**
  `getAdmeProperties` (`package.ts:133-139`): molblocks quoted without escaping internal `"` (title
  line is free text), SMILES and the column name unquoted (CXSMILES contains commas; a column named
  `smiles, canonical` corrupts the header) — one bad row fails the whole request with a pandas
  ParserError. `DG.DataFrame.toCsv` exists and this package already uses it (`admetica-utils.ts:462`).

- **AD-07 [resolved] — Container: the model checkpoint is reloaded per 1000-row batch.**
  `predict` loops models × batches (`Admetica.py:151-155`), each iteration re-running
  `load_from_checkpoint` (:82), constructing a new `pl.Trainer` (:112-118) and re-scanning the model
  dir (:128) — 5000 molecules × 1 model = 5 full model loads; the 1000-row outer batch is smaller
  than the 512 dataloader batch it wraps, so the batching buys nothing.

- **AD-08 — Table lookup by name can restyle an unrelated user table.** *(deferred 2026-09-09: skipped for now per review)*
  Views are resolved via `getTableView(name)` ("returns the FIRST TableView", js-api shell.ts:312);
  the single-molecule preview frame is named `'table'` (`constants.ts:8`, `admetica-utils.ts:387`) and
  `addColorCoding(..., true, props)` (:388) then writes `isTextColorCoded = true` onto any open user
  table that happens to be named `table` (:143-145). Unguarded null at :289 and :454.

- **AD-09 [resolved] — App init order: globals read before they're loaded.**
  `processFileData` reads `properties` at `admetica-app.ts:38` (and transitively at :47-48) before its
  own `await setProperties()` at :50. Not an unconditional crash (the Sketch-pane path initializes
  first), but a real race if the user clicks the File tab early — and the app silently inherits a
  leaked custom template (AD-02).

- **AD-10 [resolved] — Container model resolution is order-dependent.**
  `find_model` picks the first `.ckpt` from an **unsorted** `os.listdir()` substring match
  (`Admetica.py:127-133`); an empty/unknown model name matches the first file unconditionally
  (reachable from AD-03). Shipped template names are currently disjoint, so ambiguity is latent.

- **AD-11 — Preview pane writes persistent grid state and dataframe tags.** *(attempted 2026-09-10, reverted pending visual verification: candidate redesign = render via a detached `DG.Grid.create(result)` over the 1-row prediction df instead of writing settings/tags onto the real grid column)*
  Opening the context-panel Summary pane assigns pie settings to the molecule column's persisted
  `GridColumn.settings` (`admetica-utils.ts:463`) and writes `.vlaaivis-metadata` tags onto the real
  dataframe's columns (:204 via :461) — both save with the project. Crashes are caught (:436-440),
  and the settings write is load-bearing for the current renderer design — needs a redesign, not a
  hotfix.

## Low

- **AD-12 [resolved]** — First pie chart is named "piechart (1)": `pieChartIdx` starts at 1 so the
  `=== 0 ? 'piechart'` branch is dead (`admetica-utils.ts:245,249`; mirrored `admetica-form.ts:113`);
  `tablePieChartIndexMap` is keyed by table name and never cleared, so counters grow across reopens.
- **AD-13 [resolved]** — Tests: both benchmarks pass args positionally wrong (`admetica-tests.ts:129-130,
  141-142` — the string[] lands on `molecules`, TypeError before any work) — benchmark coverage is
  non-functional; 'Container. Post request' (:46) passes the not-yet-assigned module `molecules`
  (survives only because the TS body ignores `table`); `timeout: 10000000000` (≈116 days, :133).
- **AD-14 [resolved]** — `getModelsSingle` swallows errors with a commented-out `console.log`
  (`admetica-utils.ts:409-413`); `admetica-app.ts:97` renders an Error object via `ui.divText(e)`;
  commented-out init decorator left in `package.ts:22-23`.
- **AD-15 [resolved]** — Editor polish: template-save validator warns but Save still overwrites; `rgbToHex`
  drops alpha; hardcoded category colors (`constants.ts:75-81`) and `'black'` fallback
  (`admetica-utils.ts:404`). *(fixed 2026-09-10 except category colors — by design: a serialized
  chart palette for canvas rendering, not CSS styling)*

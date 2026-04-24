# Manual Test Plan (GROK-19984)

Branch: `vmakarichev/grok-19984-diff-studio-centralized-hub`. Main changes: new `hub.ts` file, new `export/export-dialog.ts` module, substantially extended `app.ts`, `package.ts`, `utils.ts`, `ui-constants.ts`, updated CSS, and a set of new icons.

## Feature 1. Hub — main landing page of Diff Studio

**Where:** `Apps > Compute > Diff Studio` launches not the editor / last model, but the Hub (`hub.ts`).

| # | Scenario | Expected result |
|---|---|---|
| 1.1 | `Apps > Compute > Diff Studio` (after `wasProcessed`, no deep-link) | Hub opens with header, description, Diff Studio icon and `Refresh / Create / Upload` buttons |
| 1.2 | Hub for a new user with no saved models | **Templates** (3 cards) and **Library** (9+ built-in) sections are visible; **Recent** is hidden or empty |
| 1.3 | `Refresh` button (sync icon) | Hub re-renders, sections rebuild, Recent picks up the latest list |
| 1.4 | `Create` button | Editor opens with the basic template |
| 1.5 | `Upload` button → pick `.ivp` from disk | File is read and opened in the editor via `runSolverApp(equations, FROM_FILE)` |
| 1.6 | `Upload` → cancel file picker | Nothing happens, no errors |
| 1.7 | Double-click on a Templates card (`Basic / Advanced / Extended`) | Corresponding template launches |
| 1.8 | Double-click on a Library card (e.g. `PK-PD`) | Model launches |
| 1.9 | Hover over a card | Tooltip with icon, model name, and description |
| 1.10 | Cursor over a card | `pointer` (not `default`) |
| 1.11 | Hub with deep-link `/apps/DiffStudio/Templates/Basic?...` | Hub is NOT shown; the model opens directly (see `isDeepLink` in `package.ts`) |

## Feature 2. Hub card context menus

**Where:** right-click on a card in any Hub section.

| # | Scenario | Expected result |
|---|---|---|
| 2.1 | Right-click on a Templates / Library (built-in) card | Menu: `Run`, `Copy link`, `Help` |
| 2.2 | `Run` → check that the model launches | Same effect as double-click |
| 2.3 | `Copy link` → paste into address bar | Opens the correct model (URL like `/apps/DiffStudio/Templates/Basic?...`) |
| 2.4 | `Help` | Datagrok side help panel opens on the relevant page |
| 2.5 | Right-click on a Recent card (built-in) | Menu: `Run`, `Copy link` (no `Help`) |
| 2.6 | Right-click on a Recent card (custom file) | Menu: `Run`, `Copy link` (URL like `/file/system.appdata/diffstudio/...`) |
| 2.7 | Right-click on a custom Library card | Menu: `Run`, `Copy link`, `Help`, `Settings...` |
| 2.8 | `Settings...` with no help URL set | `Model settings` dialog with empty `Help link` field; `Save` button is disabled |
| 2.9 | `Settings...` enter URL → Save | URL is stored in `external-models.json`, Hub re-renders, `Help` now works |
| 2.10 | `Help` for a custom model with no help URL | Warning "Help link is not defined..." |

## Feature 3. Save to Library — user models in the Library

**Where:** `layer-plus` icon on the editor's top panel (HINT `Save model to Library`).

| # | Scenario | Expected result |
|---|---|---|
| 3.1 | Open Basic template → Save to Library | Message "Saved to Library as Basic.ivp", the file appears in `System:AppData/DiffStudio/library/`, `external-models.json` is updated |
| 3.2 | Save the same model again | Name gets a `(1)`, `(2)` suffix via `unusedFileName`. The original file is not overwritten |
| 3.3 | Save a model with invalid formulas | Error via `tryToSolve` → `processError`, the file is NOT created |
| 3.4 | Save a model that was opened from the Library (`fromFile === true`) | Re-registration in the manifest without copying the file. Message `Already in Library` or `Registered in Library` |
| 3.5 | After Save to Library | The model appears on the Hub in the Library section, in Browse > Apps > Compute > Diff Studio > Library, and in Open menu > Library |
| 3.6 | Model name with unicode / special characters (`#name: My model!`) | `sanitizeModelFileName` reduces it to ASCII (falls back to `model` if nothing usable remains) |
| 3.7 | Delete the file manually via platform Files → Refresh Hub | Card disappears from the Library |

## Feature 4. Browse Tree (Templates / Library / Recent)

**Where:** `Browse > Apps > Compute > Diff Studio`.

| # | Scenario | Expected result |
|---|---|---|
| 4.1 | Expand the `Diff Studio` node | Subnodes: `Templates`, `Library`, `Recent` |
| 4.2 | Single-click on the `Templates` node | After 250 ms (`DBL_CLICK_DELAY`) a preview gallery of template cards is shown |
| 4.3 | Double-click on the `Templates` node | Preview is cancelled (`previewTimer cleared`); a regular View with the same cards opens |
| 4.4 | Expand `Library` → list of built-in + custom models | Custom models are read from `external-models.json` |
| 4.5 | Save a model to the Library from the editor | Via `LIBRARY_CHANGED_EVENT` the folder re-renders automatically; the new model appears without manual reload |
| 4.6 | Single-click on a model item in the tree | After 250 ms a preview of the model opens (without a full launch) |
| 4.7 | Double-click on a model item | Full launch via `runSolverApp` |
| 4.8 | Up/Down arrow keys on the tree | Preview updates on every navigation step; focus does NOT escape into the preview (`guardTreeFocus` 1.5 s) |
| 4.9 | Scroll the tree → click an item below the viewport | After the preview opens, scroll position is preserved and focus returns to the tree |
| 4.10 | Right-click on a built-in tree item | Menu `Run`, `Copy link`, `Help` (`Help` is absent in the Recent section) |
| 4.11 | Right-click on a custom Library tree item | Menu `Run`, `Copy link`, `Help`, `Settings...` |
| 4.12 | Open `Recent` containing an entry that points to a deleted file | The entry is skipped (`getCachedFileInfo` returns null) |

## Feature 5. Export to LaTeX / Markdown — dialog

**Where:** `file-export` icon on the editor's top panel. Implementation — `export-dialog.ts`.

| # | Scenario | Expected result |
|---|---|---|
| 5.1 | Click the Export icon with no model loaded | Message `No model loaded`, dialog does NOT open |
| 5.2 | Open a model → Export icon | Modal dialog 960×640 opens, title `Export as LaTeX`, format `latex`, preview with stex highlighting via CodeMirror |
| 5.3 | Switch Format → `markdown` | Title changes to `Export as Markdown`, preview re-generates, highlighting switches to markdown, the `Standalone` field hides, file name extension changes `.tex → .md` |
| 5.4 | Switch back → `latex` | Standalone is visible again, extension changes `.md → .tex` |
| 5.5 | Cog icon (Customize export) | Secondary options show/hide: `Title & description`, `Initial conditions`, `Parameters`, `Constants`, `Multiplication` |
| 5.6 | Turn off `Title & description` | Preview updates without the metadata section |
| 5.7 | Turn off `Initial conditions` | No initial-conditions table |
| 5.8 | Turn off `Parameters` / `Constants` | Corresponding tables disappear |
| 5.9 | Toggle `Compact` | LaTeX/Markdown without section headings, with inline initial conditions |
| 5.10 | `Multiplication` `a · b` ↔ `ab` | LaTeX toggles between `\cdot` and juxtaposition |
| 5.11 | Toggle `Standalone` (LaTeX) on | Preview adds `\documentclass{article}`, `\usepackage{amsmath, amssymb, booktabs}`, `\begin{document}…\end{document}` |
| 5.12 | `Download` button | File is downloaded with the name from the input; contents identical to the preview (including standalone wrapping) |
| 5.13 | Clear the `File name` field | `Download` button becomes disabled |
| 5.14 | Copy icon | Text is copied to the clipboard; a `Copied` badge appears next to it for ~1.5 s |
| 5.15 | Copy LaTeX → paste into Overleaf, compile with pdflatex | Document compiles without edits (with Standalone=on) |
| 5.16 | Copy Markdown → paste into Obsidian/GitHub | Math blocks render |
| 5.17 | Broken IVP (e.g. delete `#equations`) | Preview shows `// Error: <message>`, copy icon hidden, Download disabled |
| 5.18 | All `Include` toggles off | Preview is empty → copy/download disabled |
| 5.19 | Open/close the dialog several times | No memory leaks; CodeMirror is cached (`_cm`), repeated opens are faster |

## Feature 6. Open menu on the top panel

**Where:** folder icon `<i class="fas fa-folder-open">` on the left of the editor ribbon.

| # | Scenario | Expected result |
|---|---|---|
| 6.1 | Click on Open | Menu: `Import...`, `My models`, `Recent`, `Templates` (group), `Library` (group) |
| 6.2 | `Import...` → pick a `.ivp` file | Current model is overwritten with the file contents |
| 6.3 | `Library` group | Contains 9 built-in + all custom models from `external-models.json` |
| 6.4 | Click on a custom model in the `Library` group | Model loads, path is updated (`mainPath = file/...`), entry is added to Recent |

## Feature 7. .ivp file preview (file viewer)

**Where:** double-click on a `.ivp` file in `Browse > Files`.

| # | Scenario | Expected result |
|---|---|---|
| 7.1 | Open the first `.ivp` file from the file tree | Preview with full simulation opens |
| 7.2 | Open a second `.ivp` file right after the first one | Uses `file.fullPath`, not `window.location.href`; URL parameters from the first model are NOT applied to the second |
| 7.3 | Deep-link with parameters to the first file (`?param=...`) | Parameters are applied only on the initial load (`isStartingUriProcessed`) |

## Feature 8. Multistage simulation (diff-grok pipeline)

**Where:** new integration via `getIvp2WebWorker` / `getPipelineCreator` / `applyPipeline`.

| # | Scenario | Expected result |
|---|---|---|
| 8.1 | Run `Acid Production` (`#update` blocks) | Result contains an additional `Stage` column; stage separation is visible on the line chart as `segmentColumnName` |
| 8.2 | Run `PK-PD` (`#loop` for dosing) | Dose series is applied correctly |
| 8.3 | Compare numerical result against master branch (Robertson, Pollution) | Match within tolerance (best-effort: visually identical graphs) |
| 8.4 | Model with `meta.solver: {method: 'cvode'}` | CVODE is used, computation time stays within the limit |
| 8.5 | Set `maxTimeMs: 10` for a heavy model | Abort fires, a warning is shown |

## Feature 9. Caching and invalidation

| # | Scenario | Expected result |
|---|---|---|
| 9.1 | Open the Hub several times in a row | Recent table and `external-models.json` are read once (`prefetchRecentModelsTable`, `prefetchExternalLibraryEntries`) |
| 9.2 | Save to Library → Open menu | Library group is refreshed via `LIBRARY_CHANGED_EVENT` |
| 9.3 | Settings... change help URL → Hub Refresh | Custom model uses the fresh help link |

## Feature 10. Regressions (must verify)

| # | Scenario | Expected result |
|---|---|---|
| 10.1 | Run the demo `runDiffStudioDemo` | Works as before |
| 10.2 | Launch via decorator-models (`pkPdNew`, `Bioreactor`, `acidProduction`, `pollution`) | Old entry points open models via `Model.run()` (now returns the promise directly) |
| 10.3 | F5 / Refresh in the editor with modified formulas | Applies changes as before |
| 10.4 | Edit toggle ON → edit formulas → Refresh | Works; the Export icon is also available (it lives in the third ribbon panel) |
| 10.5 | Sensitivity / Fit buttons on the solve panel | Open their respective views |
| 10.6 | Save (big button) → My files | File is saved to `Home/`, Recent is updated |
| 10.7 | Download (down arrow) on the ribbon | Downloads a `.ivp` file (does NOT trigger the new Export dialog — these are different buttons) |
| 10.8 | Export to JS (`</>`) | Still creates a script |
| 10.9 | Lookup table (`#meta.inputs`) in `Bioreactor` | Works |
| 10.10 | Drag & drop a `.ivp` file into the browser | Opens via `ivpFileHandler` |

## Smoke set (minimal path for PR review)

1. `Apps > Compute > Diff Studio` → Hub is shown
2. Double-click `PK-PD` → model opens, graphs are present
3. Edit ON → change `final = 50` → Refresh → graph changes
4. Export icon → switch format, toggle Compact and Standalone, Download `.tex`, Copy → Markdown
5. Save to Library → Hub Refresh → new model visible in the Library
6. Browse > Apps > Compute > Diff Studio > Library → new model is there, right-click → Settings → set help URL
7. Hub: right-click on the custom card → Help → opens the configured URL
8. Close the app, open it again → Recent contains PK-PD and the custom model

## Additional edge cases

- Very large models (`Pollution`, 25 reactions) in the Export dialog — verify no noticeable lag on regenerate
- Behavior when there are no write permissions on `System:AppData/DiffStudio/library`
- Concurrent Save to Library (two tabs) — make sure `external-models.json` is not corrupted

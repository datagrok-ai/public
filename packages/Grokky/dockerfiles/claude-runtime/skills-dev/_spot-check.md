# Live spot-check plan — Grokky skills v1

One-line prompts per helper. Run each in the Shell AI panel (or any TVAIPanel)
on a real Datagrok session with a TableView open. Pass = the helper is reached
**and** the visible effect matches intent.

## Build & deploy

```bash
# 1. Browser bundle (grokky-api.ts changes)
cd public/packages/Grokky
npm run build && grok publish

# 2. Claude runtime container (SKILL.md changes, system prompt)
# Use whichever pipeline normally rebuilds the claude-runtime image.

# 3. Open a Datagrok session with a chemistry table:
#    `demog.csv` for non-chem prompts, any SMILES-bearing table for chem prompts.
```

## datagrok-df-and-columns

| Helper | Prompt | Pass criteria |
|---|---|---|
| `findColumns` | find the molecule column | returns hits with `column.name`, semType=Molecule scored ~1.0 |
| `describeColumn` | describe the age column | JSON with type, length, missing, numerical {min,max,avg,stdev} |
| `addColumn` | add a column 'logMW' computed from MW | new column appears in grid; type = double; values computed |
| `setColumnMeta` | set the friendly name of pIC50 to 'Potency' | grid header now reads "Potency" |
| `colorCode` via `setColumnMeta` | color-code age red-to-green from 19 to 70 | grid cells gradient red→green |
| `cloneDf` | make a copy of the filtered rows | new DataFrame returned; original untouched |
| `removeColumns` | remove the columns _tmp1 and _tmp2 | both columns gone from grid |
| `topCategories` | top 5 categories of the race column | returns array of {value,count} sorted desc |

## datagrok-filtering

| Helper | Prompt | Pass criteria |
|---|---|---|
| `filterRows` (range) | filter to rows where age > 40 | grid shows only rows with age > 40 |
| `filterRows` (substructure SMILES) | filter the molecule column by substructure CCO | grid filtered; molecule column shows matches; NOT zero matches |
| `filterByPredicate` | filter to MW < 500 AND LogP < 5 | both conditions visible in filtered grid |
| `clearFilter` | clear the filter | all rows visible again |
| `invertFilter` | invert the current filter | previously-hidden rows now visible, vice versa |
| `filteredDf` | show only the filtered rows as a new table | NEW TableView opens with subset |
| `dropRows` | remove all rows where category is X | rows count decreases; original df mutated |
| `filterSubstructure` (clear) | remove the substructure filter only | filter cleared; other filters survive |

## datagrok-selection

| Helper | Prompt | Pass criteria |
|---|---|---|
| `selectRows` (predicate) | select rows where age > 40 | selection count matches; rows highlighted in grid |
| `selectRows` (indices) | select rows 0 through 10 | exactly 11 rows selected |
| `clearSelection` | deselect everything | selection count = 0 |
| `invertSelection` | invert the current selection | selection flips |
| `selectAll` | select all rows | trueCount = rowCount |
| `selectedDf` | show only the selected rows as a new table | new TableView with subset |
| `filterFromSelection` | filter to just the selected rows | filter mask = selection; non-selected hidden |
| `selectionFromFilter` | select all currently-visible rows | selection.trueCount = filter.trueCount |
| `setCurrentRow` | go to row 5 | current row indicator moves to row 5 |
| `describeSelection` | what's selected? | returns {count,total,indexes,currentRowIdx,sample} |

## datagrok-viewers

| Helper | Prompt | Pass criteria |
|---|---|---|
| `addViewer` (scatter) | add a scatter plot of height vs weight | scatter plot docks in view; correct axes |
| `addViewer` (histogram) | add a histogram of age | histogram appears |
| `addViewer` (color by) | scatter MW vs LogP, color by activity | colorColumnName set |
| `configureViewer` | turn on the regression line on the scatter plot | line appears; viewer not re-created |
| `findViewer` | find the scatter plot | returns viewer handle; skipping the grid |
| `findViewers` | list all charts | returns array of non-grid viewers |
| `closeViewer` (type) | close the scatter plot | scatter plot removed; histogram/others remain |
| `closeAllViewers` | close all viewers except the grid | only grid remains; filter panel may remain |

## datagrok-grid-customization

| Helper | Prompt | Pass criteria |
|---|---|---|
| `applySort` (single) | sort by activity descending | grid re-orders by activity desc |
| `applySort` (multi) | sort by class asc then activity desc | both columns drive sort |
| `clearSort` | clear the sort | grid back to natural order |
| `configureGrid` (visibility) | show only SMILES, MW, LogP, activity | grid shows exactly those 4 |
| `configureGrid` (widths) | make the molecule column 300 pixels | column width = 300 |
| `configureGrid` (formats) | format IC50 to 2 decimals | IC50 cells render as "0.00" |
| `pinColumn` | pin SMILES to the left | SMILES becomes leftmost, sticky on scroll |
| `unpinColumn` | unpin SMILES | SMILES un-pins |
| `colorCode` (linear) | color activity from red to green | gradient renders in grid |
| `colorCode` (off) | turn off color coding on activity | gradient removed |
| `resetGrid` | reset all grid customization | visibility/widths/sort/colors all default |

## Cross-skill recipe (the multi-step demo from the RFC)

| Prompt | Pass criteria |
|---|---|
| Pin SMILES on the left, hide the index column, color activity red-to-green, sort by activity desc | all four effects visible together in one go |

## What "fail" looks like for each

- Helper not reached → Claude wrote raw js-api code. Check the `datagrok-exec` block in chat.
- Helper reached, no visible effect → likely a runtime exception. Open browser DevTools console, check for "Uncaught" or stack traces. Capture and report.
- Helper reached, partial effect → behavioral bug. Note which dimension fails (e.g., "applies sort but loses width customization").

## Reporting back

For each failure, capture:
1. The prompt sent
2. The exec block Claude emitted (visible in chat)
3. The console output (DevTools)
4. The visible grid/viewer state (screenshot if useful)

Send these back and I'll diagnose. Most likely categories:
- **Skill defect** — Claude reaches a different helper than expected
- **Helper bug** — the helper itself crashes or behaves wrong
- **API drift** — `datagrok-api` shape changed since research

## Per-turn view context refresh (Phase 5 side-effect)

You should see Claude reference current state (selection count, filter count,
current row, semType-tagged columns) **in every turn**, not just turn 1. If
Claude says "I don't know what's selected" mid-conversation, the per-turn
prepend isn't wired. Look at the prompt that hits the runtime (visible in
`claude-runtime` logs) — first line should be `Table "..." (...) ...\nState: ...`.

# Legend — what the features cannot say yet

Written during the Grid / Filter panel / Legend gap round of 23 Oct 2026, against the old
`playwright-public/legend` specs, `files/TestTrack/Viewers/Legend/*.md`, the legend's DOM contract
(`core/client/d4/lib/src/viewer_base/legend_mixin.dart`, `widgets/legend/`) and the library's legend
steps (`bindings/tiers/viewers/steps.ts`, `src/runtime/viewer-legend.ts`); trimmed in the BDD PR round of 8 Nov 2026.

## Resolved

- Down and Up move a sole selection to the next and the previous category again (public#2548, core#1159: the handler
  the legend refactor dropped is back in `legend_interaction.dart`): `legend-interaction.feature`.

## Still open

- **Per-section readings**: a scatter plot coloured by one column and marked by another has two legend sections; the
  legend publishes one total (`data-legend-items`); needs `data-legend-items-<section>` on the legend root.
- **The inner splitter between two sections** (`LegendView.initInnerSplitter`): needs a name (`legend-inner-splitter`)
  and per-section heights.
- **The markers column selector in the legend header**: re-picking DIS_POP lists 9 items, not the 10 categories (4 RACE
  + 6 DIS_POP), while a later re-layout lists 10 (measured 2026-10-24); needs a ruling first.
- **A label that does not fit**: whether an item's label is cut; needs a `data-legend-truncated` count or a per-item
  overflow flag.
- **Scrolling a long list**: the list renders a window of its items; needs the first rendered item key
  (`data-legend-first-item`) and a step that scrolls the legend list by a page.

## Won't do

- Up on its own from Black (public#2548): the Down-and-Up scenario claims it.
- The colour column renamed or removed under a selection: the shared viewer rule (ApiTests `AI: Viewers: Shared
  rules`, Dart `Column removal`); that a rename releases the selection is incidental.
- A filter that removes every row under a selection: the category-dropping filter in "A selection stays through a row
  filter…" claims the same; no rows is the shared viewer rule.
- A 300-character label, 20001 rows with Position Auto, a one-value colour column: cosmetic or not UI checks.

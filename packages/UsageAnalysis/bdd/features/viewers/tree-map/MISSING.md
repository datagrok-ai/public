# Tree map — what the features cannot say yet

Written during the heat map / 3D scatter plot / tree map gap round of 19 Oct 2026, against
the viewer's "Automation surface" (`core/client/d4/lib/src/viewers/tree_map/CLAUDE.md`) and the library's
rules (`public/libraries/bdd/CLAUDE.md`). `tree-map-properties.feature` holds what could be said:
Show Selected Rows, Row Source Selected, a bar chart filtering the map (public#935), Auto Layout and the
selectors (public#2131), and the categorical gate.

## Still open

- **Picking a column in an on-viewer split selector.** The split selectors are plain `<select>` elements
  (`class="d4-column-selector-aggregation d4-column-selector-tree-map"`, `TreeMapCore._addSplit`), not
  `ColumnComboBox`es: neither `user picks {string} in the {string} column selector of {widget}` nor
  `user picks {string} in the column selector at the {string} area of {widget}` reaches them, because both
  look for `.d4-column-selector` (`gestures.ts` `select`, `widgets.ts` `pickInAreaSelector`), and the
  selects carry no name a `{element}` could resolve. Blocked, so not written:
  - picking a column in the trailing empty selector (`split selector 2` on a one-level map) adds a level:
    `levels` 1 → 2, `selectors` 2 → 3, `split columns` "RACE, SEX";
  - setting a middle selector to empty drops it and every level after it (`_splitChange` removes the
    selectors past the emptied one): `levels` 2 → 1 on "RACE, SEX".
  Needs either a name on each select (`split selector <n>`, the area's own name) that `user selects
  {string} in {element}` can resolve, or the area-based pick learning a native `<select>` under the area.
  The features set `splitByColumnNames` as a property instead, which is the other half of the same round
  trip (`onLookChanged` rebuilds the selectors).
- **`{widget} should report the error {string}` on the tree map.** The generic step reads
  `WidgetStatus.error`; the tree map's status leaves it null and reports the gate as the `error`
  *reading* (`values['error']`, `_TreeMapStatusFeature._error`), so the step read "" on a table with no
  categorical column. The feature claims `the "error" reading of tree map viewer should be "Requires
  categorical columns"`. Either the status fills `error:` as well, or the step falls back to the reading.
- **The selectors' absence comes from the same getter the reading does.** `tree_map_status.dart:42`
  (`core/client/d4/lib/src/viewers/tree_map/features/`) reports `split selector <n>` and `color selector`
  only `if (plot.showColumnSelectionPanel)` — the computed getter (Auto Layout and the 200 by 200 rule), not
  the DOM. So the Auto Layout scenario's "should not have a `split selector 1` / `color selector` area"
  claims restate the `selection panel shown` reading and would pass with the selectors still on screen
  (a `htmlSetDisplay` that did not run). Fix: drop the `if`, so each area is reported only while
  `_controlRect`'s `htmlIsElementVisible` finds the control shown.
- **No selection highlight cannot be claimed on the RACE colours.** `tree map viewer should show no
  selection highlight` reads 275 587 px in the selection hue on demog-1000 split by RACE with Row Source
  Selected and no band drawn (measured 2026-10-19): the category colours themselves fall in the selection
  hue, so the claim does not tell a band from no band and is not written. Row Source Selected's "no band"
  is claimed by the absent `selection of <path>` areas only. Needs either the band counted in the status
  (a `bands drawn` reading) or a pixel claim restricted to the band's colour as drawn.

## Answered

- **What a person sees on a table with no categorical column** (lead, 2026-10-19): intended. Add Viewer
  does not offer the tree map for such a table — the item is disabled with the tooltip "Requires
  categorical columns" — so a person never reaches the blank map. Only a viewer added through the API
  draws it, with the gate in the `error` reading, which is what the feature claims.

## Won't do

- **The selection band's height** (`selection of <path>` as the leaf's selected share): the area is the
  renderer's own condition and rectangle; the features claim it present or absent and pair the toggle
  with a repaint, not a pixel scan of the band.
- **What a layout stores for Show Selected Rows / Auto Layout**: a layout round trip is ApiTests',
  `tree-map-color-size-selection.feature` already claims the one the md asks for.

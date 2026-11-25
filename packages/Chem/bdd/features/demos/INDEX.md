# Demo clips

Chem's showcase clips, in the order their reel cuts them. Each clip is a scenario of a feature in this folder and a
test of what it shows. `MSYS_NO_PATHCONV=1 npx grok-bdd guide features/demos/<feature> --gif` in `bdd/` films it into
`guides/<feature slug>/<scenario slug>/`; `py public/libraries/bdd/tool/guide-render.py --reel <out.mp4> <clip dirs…>`
cuts the clips together. The `/bdd-answer` skill, section 5, is the recipe.

## Reactions drawn by Crux (2026-11-13)

Chem's reactions in Crux mode, over the showcase table `fixtures/BDD-reactions-showcase.csv` (named reactions with
their agents and products, four routes with a ring scaffold, seven of Chem's enumeration templates).

| # | Clip | Feature |
|---|---|---|
| 1 | Reactions drawn like the molecules beside them | `reactions-crux.feature` |
| 2 | Routes drawn as whole chains | `reactions-crux.feature` |
| 3 | The Reaction pane maps a route step by step | `reactions-crux.feature` |
| 4 | Templates drawn as generic schemes | `reactions-crux.feature` |
| 5 | Switching the Chemical Engine | `reactions-crux.feature` |

## Reactions edited in Crux Sketch (2026-11-14)

The Reactions app's sketcher (land 5; the app is the Reaction Sketcher since 2026-11-16): the text and the drawing
following each other.

| # | Clip | Feature |
|---|---|---|
| 6 | The Reactions app draws what you type, and writes what you draw | `reactions-edited-crux.feature` |

Retired 2026-11-16 for the clips below: "A reaction edited in its cell" (the editor opened from a cell, edited and
written back, is clip 15's) and "A route the editor does not open says why" (the editor opens every route since
crux-core 15b1bbd65; Chem's package tests keep the refusal).

## The Reaction Sketcher and reaction tables (2026-11-16)

The reaction work of the run of 2026-11-14 to 16, as a chemist meets it: the Reaction Sketcher (Apps > Chem >
Reactions) and what its canvas does, mapping, Export for a journal and Add to table; then reaction cells, the reaction
editor opened from a cell and ChemDraw pages. Each Sketcher clip opens the app on a reaction of its own; the tables are
Chem's `files/tests/reactions-crux.csv`, `files/reactions-cdxml.csv` and the showcase table, opened from the project.
Its reel is cut from these ten alone.

| # | Clip | Feature |
|---|---|---|
| 7 | The Reaction Sketcher, on its example route | `reaction-sketcher.feature` |
| 8 | Building a route on the canvas | `reaction-sketcher.feature` |
| 9 | Conditions over an arrow | `reaction-sketcher.feature` |
| 10 | A catalytic cycle | `reaction-sketcher.feature` |
| 11 | Map numbers through a route | `reaction-sketcher.feature` |
| 12 | Export for a journal | `reaction-sketcher.feature` |
| 13 | The route into a table | `reaction-sketcher.feature` |
| 14 | Reaction cells, from SMILES and from ChemDraw | `reaction-tables.feature` |
| 15 | The reaction editor, from a cell | `reaction-tables.feature` |
| 16 | A ChemDraw page, through the editor and back | `reaction-tables.feature` |

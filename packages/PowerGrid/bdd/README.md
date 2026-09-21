# Power Grid behavioral tests

Gherkin features under `features/`, compiled by `@datagrok-libraries/bdd` (`public/libraries/bdd`)
into the Playwright specs under `generated/` — committed, never edited by hand. They replace the
TestTrack Forms viewer specs (`packages/UsageAnalysis/files/TestTrack/Viewers/FormsViewer/`), whose
assertions each feature carries or strengthens (the review record is the `/bdd-translate` skill
under `public/.claude/skills/`). Five features, 30 scenarios, 25 s in one run on one worker
(2026-09-09).

The subject is the **Forms viewer** — the PowerGrid viewer that shows several records side by side
(`libraries/utils/src/viewers/forms-viewer.ts`, registered as `Forms` in `src/package.ts`). It is
not the core Form viewer, whose features live in `packages/UsageAnalysis/bdd/features/viewers/form/`.

| Feature | Scenarios | What it claims |
|---|---|---|
| `forms-core` | 9 | one card for the current row and one per selected row that passes the filter, in table order; the grid's sort mirrored onto the card order and the header's arrow; Sort By overriding it without touching the grid's own sort; Use Grid Sort off stopping the mirror (GROK-20380, fixed here); the sort cycle a double-click on any label advances; Pin Row moving a card into the pinned pane, Unpin Row bringing it back, and the warning a non-unique value earns |
| `forms-interactions` | 10 | a click on a card makes its row current, a click on a field sets the current cell, a click on a header label sets the current column, Control toggles a row's selection, Shift selects every row up to the card's and Control+Shift clears them, hovering a card or a grid row fills the mouse-over card and leaving empties it, and each show toggle drops its card |
| `forms-fields` | 6 | the Fields property picks the columns and their order, the header's remove icon drops one, a column renamed to a `~` name leaves the set, an empty set draws no card and says nothing, a column removed from the table takes its field with it, and a named number format reaches the float fields only |
| `forms-persistence` | 2 | a non-default field set, the Sort By column and a row pinned by value survive a layout saved to the server and re-applied over a view whose Forms viewer was closed and a histogram added, and a project round-trip |
| `forms-presentation` | 3 | a text column is an input field; a column's colour coding paints the field's background on the current card and on the selected rows' cards, Color Code gates it and removing the coloring clears it |

From a fresh checkout of `public`, against a local stand on `http://localhost:8888` with this
package published:

```bash
cd public/libraries/bdd && npm ci && npm run build   # the library (a path dependency; dist/ is not committed)
npx playwright install chromium                      # its browser, once per machine
cd ../../packages/PowerGrid && npm install           # the package
npx grok-bdd link                                    # ONE Playwright: the library's copy into node_modules (redo after every npm install)
npx webpack && grok publish localhost                # the features read the viewer's own status
npx grok-bdd run --reporter=list                     # compile --check, then Playwright (one page per worker; --workers N)
npx grok-bdd run generated/viewers/forms             # one folder
```

The stand needs a platform from `core` at or after 2026-09-08 (the grid's `getWidgetStatus`, from
which the features read `header <column>`, `cell <row> of <column>`, `sort column` and `current
column`) and demog-1000 under `System:DemoFiles` — every literal in the features is a value of that
file. The viewer's own status is what the features read: `cards`, `records shown`, `pinned
records`, `fields`, `fields shown`, `sort column`, `sort direction`, `current record`,
`mouse-over record`, `<COLUMN> of card <k>` (and of `current card`, `mouse-over card`, `pinned card
<k>`), `field kind of <COLUMN>`, `width`/`height`/`background of <COLUMN> of card <k>`; the hit
areas are `card <k>`, `current card`, `mouse-over card`, `pinned card <k>`, `field <COLUMN> of card
<k>`, `label <COLUMN>`, `remove <COLUMN>` and `sort indicator <COLUMN>`. It is built in
`libraries/utils/src/viewers/forms-viewer-status.ts`, so publishing a change to it means building
`libraries/utils` before this package.

The features change the selection, the filter, the grid's sort, the field set and the pinned rows,
and put each back; the persistence feature's project is deleted at the feature's end. The molecule
and curve scenarios of the old specs (the renderer-size ladder, the substructure filter) are not
translated: they need Chem and Curves published, and this stand has neither.

Editing: change a feature, `npx grok-bdd compile`, commit the regenerated spec with it;
`npx grok-bdd list-steps` prints every phrase this package can use, including its own
(`bindings/forms.ts`: the pinned pane's visibility, a column rename, the current column, and a
click with two modifiers held).

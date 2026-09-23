# PowerPack behavioral tests

Gherkin features under `features/`, compiled by `@datagrok-libraries/bdd` (`public/libraries/bdd`)
into the Playwright specs under `generated/` — committed, never edited by hand. They translate the
TestTrack PowerPack section (`packages/UsageAnalysis/files/TestTrack/PowerPack/`), one feature per
dialog session or data source, each saying in its description what it claims and what it leaves out.

| Folder            | Features | What it claims |
|-------------------|----------|----------------|
| `add-new-column/` | add new column, formula editor, functions panel, formula refreshing, persistence over a file in My files, persistence over Northwind query results | the Add New Column dialog from the toolbar and from Edit: tooltips, resize, a formula built from autocomplete and columns dragged from the grid's header and from the dialog's column list, the input history; autocomplete on a letter, Ctrl+Space and "$", the signature on hover, column references highlighted in --blue-2 (pasted, autocompleted, dropped, and the GROK-17004 formula); functions inserted by the plus icon and by a drag, the picked column passed when its type fits, the list reordered by the picked column's type and kept by "By name"; a chain of calculated columns recalculated from the Formula pane of the context panel; calculated columns following a rename and an edit of their source and coming back from a project saved with Data sync (GROK-17109) |
| `annotation-regions/` | annotation regions, drawing, interaction, interaction on every viewer, titles, persistence, the Tools menu, the other 2D viewers | regions drawn and edited on the scatter plot, line chart and the other two-dimensional viewers, hovered and clicked, their titles and title strips, kept with the view — through the regions' hit areas and readings |
| `formula-lines/` | the Formula Lines dialog, the axis Annotations group, regressions, GitHub regressions | lines and bands added, edited, hidden and removed in the dialog and the look it writes; the Annotations group on a numeric axis; the fixed regressions, and #2487 (a datetime axis), #671 (the preview follows the viewer's axes and the selected line), #2747 (one dataframe line on a scatter plot and a line chart) |
| `home/` | home widgets, home as sharing user | the Home widgets, hiding them, Customize, a reload, the Spotlight tabs, the Community links, Usage and Reports; a second account without Usage and Reports, an unread notification on the server |
| `search/` | power search enter | Enter on the md's queries once the results have answered (`aria-busy=false`), with and without the suggestion menu, and the suggestions walked with the arrows |
| `navigation/` | direct link loading | a project opened by its direct link on a fresh page load and the same project from the Browse tree |
| `io/` | xlsx open, xlsx shared with me | a three-sheet workbook through My files, a drop, the file chooser, File > Open and Shared with me — exactly three table views |
| `enrichment/` | data enrichment | creating, applying, editing and deleting enrichments, several at once, a layout and a reopened project, the enrichments of a referenced column offered through a foreign key |

`annotation-regions/` and `formula-lines/` came from `packages/UsageAnalysis/bdd/features/viewers/`:
the feature belongs to PowerPack (the Formula Lines dialog is `PowerPack:formulaLinesDialog`).

`add-new-column/persistence-northwind` is `@full-stand`: it needs the Postgres Northwind database
(the NorthwindTest connection of dev); a stand without it runs with `--grep-invert @full-stand`.

What the stand needs: PowerPack published from the same checkout, Chem (the SPGI features pick
Molecule functions), the demo files under `System:DemoFiles/` and a "My files" share for the
current user; `home/` and `io/xlsx-shared-with-me` also the library's sharing account, `home/` an
administrator account of its own (created once per stand), `enrichment/` the `System:Datagrok`
connection.

Run from the package's `bdd` directory. PowerPack is a package of the `public/` pnpm workspace, so
it resolves the library with nothing to link:

```bash
cd public && grok setup                               # once per checkout: the pnpm workspace
cd libraries/bdd && npm run build                     # the library; dist/ is not committed
cd ../../packages/PowerPack/bdd
npx grok-bdd run generated --workers=2 --reporter=line --grep-invert @full-stand
npx grok-bdd run generated/add-new-column             # one folder
```

The bindings are the package's own screen parts and checks: `bindings/add-new-column.ts` — the
Add New Column dialog's editor, hint, column list, functions list, preview and input history; the
column list's rows found by the grid's own readings; the functions list's order and what its top
rows take; the highlighted column references and their computed color; relations between calculated
columns checked row by row; a file put into the user's My files share for the feature.

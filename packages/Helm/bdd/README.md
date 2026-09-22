# Helm behavioral tests

Gherkin features under `features/`, compiled by `@datagrok-libraries/bdd` (`public/libraries/bdd`)
into the Playwright specs under `generated/` — committed, never edited by hand. They replace the
TestTrack Helm cases (`packages/UsageAnalysis/files/TestTrack/Helm/`) and the hand-written specs
under `playwright/`, whose UI assertions each feature carries or strengthens (the review record is
the `/bdd-translate` skill under `public/.claude/skills/`). Six features, 23 scenarios, about 30 s
on three workers against dev (2026-09-22, green three runs in a row).

| Folder     | Features             | What it claims                                                                 |
|------------|----------------------|--------------------------------------------------------------------------------|
| `render/`  | renderer             | the Bio detector tags the HELM column, the grid reports `helm` as its cell type and paints monomers in their colors, before and after a scroll |
| `editor/`  | open, notation, palette | a double-click and Current Value > Edit Helm... open the editor on the cell's sequence (exact notation, monomer count, toolbar, palette, tabs); an invalid notation shows the parse error and keeps the drawing, a valid one redraws without touching the cell; undo/redo, Clean layout, the Properties tab; OK writes the exact edited notation, Cancel nothing; the palette's search, a tile arming "Next add" and a canvas click placing it, the RNA triplet builder, the empty Favorites |
| `panels/`  | properties           | the context panel's Properties pane: formula, weight and extinction coefficient of the current cell, following the current cell, and the "Too long sequence" guard over 1000 characters |
| `service/` | surface              | `Helm:getHelmHelper` exposes the methods other packages call; `Helm:getMolfiles` returns one hwe pseudo-molfile per row |

One scenario is `@known-failure`, GROK-20962: Edit Helm... opens the current row, not the cell it was picked on (`openEditor` in `src/package.ts`
reads `df.currentRowIdx`), so a right-click on another cell edits — and on OK overwrites — the
current one.

What is not here, and why (each feature says it too): the monomer tooltip on hover, since the
renderer publishes no hit area for the monomers it draws (a `grid.addStatusProvider` in
`HelmGridCellRendererBack` would give one); the helper's computations (`parse`, `removeGaps`,
`getHoveredAtom`, `createHelmInput`, `createHelmWebEditor`, the monomer-function override,
`buildMonomersFuncsFromLib`), which have no user-visible effect and are the package's own tests in
`src/tests` (`helm-helper-surface-tests.ts` covers the ones no test had); starring a Favorites tile,
which would leave a user setting behind; the extinction coefficient of the editor's Properties tab,
which shows 0 where the context panel shows 0.06 for the same sequence — the editor does not show
the full number (a display format, confirmed by the Helm owners), so the pane's value is claimed.

Run from the package directory against a stand with Helm, Bio and Chem published. Helm is a
package of the `public/` pnpm workspace, so it resolves the library with nothing to link:

```bash
cd public && grok setup                               # once per checkout: the pnpm workspace
cd libraries/bdd && npm run build                     # the library; dist/ is not committed
cd ../../packages/Helm
DATAGROK_URL=https://dev.datagrok.ai DATAGROK_SERVER=dev npx grok-bdd run --workers=3 --reporter=list
npx grok-bdd run generated/editor                     # one folder
```

The datasets are the package's own samples under `System:AppData/Helm/samples/`
(`bindings/elements.ts`), published with it; the features change their in-memory copy only and
leave nothing on the server. The editor's parts carry the hwe library's `data-testid`s, not the
u2 contract, so `bindings/elements.ts` names them, each scoped to `HELM editor`; the editor takes
Control+A for its own select-all, which is why the notation is typed through `user replaces the
text of {element} with {string}` (`bindings/steps.ts`) rather than the library's typing.

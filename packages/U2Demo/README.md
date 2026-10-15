# U2 Demo

Demo and consumer package for [`@datagrok-libraries/u2`](../../libraries/u2), the next-gen
Datagrok UI library (design tokens, headless signal-based components, platform bridge).

Open **Apps | Dev | U2 Demo**. A nav tree on the left holds 30 sub-demos in eight areas; the same
tree hangs under the **U2 Demo** node in Browse, and every leaf round-trips to the URL
(`/apps/U2demo/U2Demo/<area>/<sub-demo>`; the last segment may also be the sub-demo's label,
`/Sections%20%26%20wizard`). With nothing remembered and no path in the URL the app opens on
**Start / Overview**.

- **Start** — *Overview* (what u2 is, what each area covers, and the three sub-demos to read first)
- **Inputs** — *All inputs* (every input type, Dart `ui.input` and u2 `propertyForm` side by side),
  *Basic inputs* (signals, binding and validation), *Range slider* (`RangeSlider`),
  *Multi-select* (`MultiSelect`, `ButtonGroup`), *Async* (Combobox over local and async items,
  `AsyncView` loading states)
- **Containers** — *Layout* (splitter, accordion, breadcrumbs, panel-local toolbar),
  *Popups* (menu, context menu, tooltip, dialog)
- **Collections** — *Lists* (`VirtualList` over 100,000 items), *Trees* (`VirtualTree` with lazy
  branches and `expandPath`)
- **Display** — *Cards* (`Card` surfaces and `StatCard` KPIs), *Feedback* (`ProgressBar`,
  notification balloons, the guided `Tour`), *Tables* (`BasicTable` and the virtualized
  `VirtualGrid`), *Sections & wizard* (`Section`, `Wizard`), *Message input* (`MessageInput`)
- **Forms** — *Form*, *Property grid*, *Object form*, *Functions* (forms built from `Func`
  metadata), *FuncCalls* (the Dart-vs-u2 bench over real `FuncCall`s), *Run history* (Run and the
  saved runs of a `funcForm`)
- **Platform** — *Dataframes* (table/column/columns pickers over live frames), *Files*, *Entities*,
  *Spaces*, *Molecules* (needs the Chem package), *Bridge* (the `u2/dg` layer: `asDartInput()`,
  `leakReport()`), *Filters* (`FilterBuilder` and `FilterQueryInput` over a DataFrame, a domain
  table and entity metadata)
- **Automation** — *MSA workbench* (a small app the bdd library's own features drive), *Catalog*
  (every registered u2 component in each state its props reach, generated from the registry)

The ribbon carries the demo commands — show this sub-demo's source in the context panel, the
Inspect toggle (click any control to see its live properties there), rebuild, and a Demo tools
drop-down. Where you are and the last action show in the shell status bar.

The whole demo is one `Component.build(...)` tree: closing the view runs the platform kill-walk,
which disposes every effect and listener — verify with the leak detector on *Bridge*.

## Behavioral tests

`bdd/features/demo/**` covers every sub-demo with [`@datagrok-libraries/bdd`](../../libraries/bdd)
features; they need this package published on the stand (`npx grok-bdd run --reporter=list` from
`bdd/`). `user opens the "<label>" demo page` opens a sub-demo by its nav label and enters the
`U2 Demo` context, where phrases resolve inside the sub-demo first (`demo navigation` is the tree).
Everything else is the library's vocabulary, by what the page shows: labels, placeholders
(`"Search users…" typeahead`), and the `name = value` status lines as readouts (`value of rows
readout`; `value of first sync readout` where a name repeats).

`bdd/features/demo/catalog.feature` is generated: `node bdd/tools/gen-catalog.mjs` writes it from
`libraries/u2/manifest.json` and `src/pages/catalog-states.json` (one source with the Catalog page),
and `--check` reports drift; after a u2 manifest change, regenerate and `grok-bdd compile`. A claim a
u2 defect fails goes into `bdd/tools/catalog-known.json` (`{"<tag>": {"<step text>": "<owner: reason>"}}`),
which the generator turns into a pinned `@known-failure` and a `TODO:` line; the entry goes when u2
is fixed. `bdd/bdd.config.json` runs the layout, console-error and leak audits in report mode: the
findings are warnings and `bdd-audit` annotations in `bdd/test-results/report.json`, summed up by
`node ../../libraries/bdd/tool/audit-summary.mjs bdd/test-results/report.json`.

## Development

```bash
npm install
npm run build          # grok api && grok check --soft && webpack
grok publish dev       # publish to the dev server in debug mode
```

The u2 library is consumed as a local file dependency (`../../libraries/u2`); rebuild it with
`npm run build-all` when the library changes.

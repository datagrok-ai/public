# U2 Demo

Demo and consumer package for [`@datagrok-libraries/u2`](../../libraries/u2), the next-gen
Datagrok UI library (design tokens, headless signal-based components, platform bridge).

Open **Apps | Dev | U2 Demo**. A nav tree on the left holds 26 sub-demos in seven areas; the same
tree hangs under the **U2 Demo** node in Browse, and every leaf round-trips to the URL
(`/apps/U2demo/U2Demo/<area>/<sub-demo>`). With nothing remembered and no path in the URL the app
opens on **Start / Overview**.

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
  metadata), *FuncCalls* (the Dart-vs-u2 bench over real `FuncCall`s)
- **Platform** — *Dataframes* (table/column/columns pickers over live frames), *Files*, *Entities*,
  *Spaces*, *Molecules* (needs the Chem package), *Bridge* (the `u2/dg` layer: `asDartInput()`,
  `leakReport()`)

The ribbon carries the demo commands — show this sub-demo's source in the context panel, the
Inspect toggle (click any control to see its live properties there), rebuild, and a Demo tools
drop-down. Where you are and the last action show in the shell status bar.

The whole demo is one `Component.build(...)` tree: closing the view runs the platform kill-walk,
which disposes every effect and listener — verify with the leak detector on *Bridge*.

## Development

```bash
npm install
npm run build          # grok api && grok check --soft && webpack
grok publish dev       # publish to the dev server in debug mode
```

The u2 library is consumed as a local file dependency (`../../libraries/u2`); rebuild it with
`npm run build-all` when the library changes.

# U2 Demo

Demo and consumer package for [`@datagrok-libraries/u2`](../../libraries/u2), the next-gen
Datagrok UI library (design tokens, headless signal-based components, platform bridge).

Open **Apps | Dev | U2 Demo** for a tabbed tour:

- **Inputs** — text, bool, number, choice, multi-choice; validation; signal binding and computed readouts
- **Combobox & async** — local and async item sources over the shared `AsyncSource` contract; `AsyncView`
- **Lists & trees** — `VirtualList` with 100,000 rows, `VirtualTree` with lazy children and `expandPath`
- **Layout** — `Splitter`, `Accordion` (lazy panes), `Breadcrumbs`
- **Popups** — anchored and context `Menu`, modal `Dialog`, shared `Tooltip`
- **Form & properties** — `Form` with aggregated validity, `PropertyGrid`
- **Platform bridge** — the `u2/dg` layer: `host()` (DG.Widget lifecycle), `asDartInput()`
  (u2 input in a `DG.Dialog`), table/column pickers, `leakReport()`

The whole demo is one `Component.build(...)` tree: closing the view runs the platform kill-walk,
which disposes every effect and listener — verify with the leak detector on the last tab.

## Development

```bash
npm install
npm run build          # grok api && grok check --soft && webpack
grok publish dev       # publish to the dev server in debug mode
```

The u2 library is consumed as a local file dependency (`../../libraries/u2`); rebuild it with
`npm run build-all` when the library changes.

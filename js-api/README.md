# datagrok-api

TypeScript/JavaScript API of the [Datagrok](https://datagrok.ai) platform. Plugins (packages) import it; the
platform ships the matching runtime, so the package is an `external` in plugin bundles, never bundled.

```ts
import * as grok from 'datagrok-api/grok';   // platform services: shell, data, dapi, functions, events, ...
import * as ui from 'datagrok-api/ui';       // DOM builders, inputs, dialogs, menus
import * as DG from 'datagrok-api/dg';       // classes, enums, and types: DataFrame, Column, Viewer, ...
```

Reference: [datagrok.ai/api/js](https://datagrok.ai/api/js) (index) and the per-class pages under
[datagrok.ai/js-api](https://datagrok.ai/js-api/dg/classes/DataFrame). Runnable samples:
[public.datagrok.ai/js/samples](https://public.datagrok.ai/js/samples) (sources in
[`packages/ApiSamples/scripts`](https://github.com/datagrok-ai/public/tree/master/packages/ApiSamples/scripts)).
Developer guide: [datagrok.ai/help/develop](https://datagrok.ai/help/develop/packages/js-api).

## Where things live

| Looking for... | Start with |
|---|---|
| Tables, columns, rows, selection, filter | `DG.DataFrame`, `DG.Column`, `DG.BitSet` (`src/dataframe/`) |
| Opening tables and views in the workspace | `grok.shell` (`src/shell.ts`) |
| Calling and registering functions, scripts, queries | `grok.functions`, `DG.Func`, `DG.FuncCall` (`src/functions.ts`, `src/entities/func.ts`) |
| Server entities: users, groups, files, projects, connections | `grok.dapi.*` (`src/dapi.ts`) |
| Viewers | `DG.Viewer.scatterPlot(df, options)` and the other statics (`src/viewer.ts`); the grid in `src/grid.ts` |
| UI: layout, inputs, dialogs, menus, tooltips | `ui.div`, `ui.input.*`, `ui.dialog`, `DG.Menu`, `ui.tooltip` (`ui.ts`, `src/widgets/`) |
| Platform events | `grok.events.on*`, per-object `on*` getters (`src/events.ts`) |
| Constants and enums | `DG.TYPE`, `DG.COLUMN_TYPE`, `DG.VIEWER`, `DG.SEMTYPE`, `DG.TAGS` (`src/const.ts`) |
| Custom viewers, filters, cell renderers, object handlers | `DG.JsViewer`, `DG.Filter`, `DG.GridCellRenderer`, `DG.ObjectHandler` |
| Domain tables (entity-mapped plugin schemas) | `grok.dapi.domains.table('<schema>.<table>')` (`src/domains.ts`) |

## Three things worth knowing before writing code

- `grok.functions.call(name, params)` resolves to the single output value when the function declares one output,
  and to an object keyed by output name when it declares several.
- `grok.dapi.<name>` creates a fresh data source on every access, and `filter/order/page/by` mutate the source they
  are called on — chain them on one instance and finish with `list()`.
- Every wrapper class holds its Dart handle in `.dart`; pass the wrapper itself to API methods (they unwrap), and
  never hand a raw handle back to your own code without `DG.toJs()`.

## Building

```bash
npm install
npm run build        # tsc + webpack + bundle smoke test; also writes the browser bundle used by the platform client
npm run build-ts     # TypeScript only
```

See `CLAUDE.md` in this folder for the interop conventions, the module map, and the regeneration rules for
`src/api/*.g.ts`. Changes go in `CHANGELOG.md` under `## v.next`.

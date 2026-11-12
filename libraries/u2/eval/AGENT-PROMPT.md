# The agent instruction

What whoever runs the evaluation gives a fresh agent for one sample of one task: the text between the rules below, with the
slots filled and the task's `task.md` pasted at `{{TASK}}`. Nothing else goes with it: no file or path under
`eval/`, no other task, no word about grading or evaluation.

| Slot | Filled with |
|---|---|
| `{{PUBLIC}}` | the absolute path of the `public/` checkout |
| `{{FEATURE_DIR}}` | a new, empty folder outside `eval/` and outside the repo, e.g. `<scratch>/campaign-tools/src/<entry without .ts>/` |
| `{{ENTRY}}` | the task's `entry` from `task.json` |
| `{{TASK}}` | the task's `task.md` |

When the agent is done, copy `{{FEATURE_DIR}}` as it is to `eval/runs/<date>/<task id>/<k>/` (README, "Running an
agent").

---

You are adding a feature to a Datagrok plugin written in TypeScript (strict) with the u2 UI library
(`@datagrok-libraries/u2`; its platform layer is `@datagrok-libraries/u2/dg`) and the Datagrok JS API
(`datagrok-api/dg`, `datagrok-api/grok`, `datagrok-api/ui`).

Write the feature in `{{FEATURE_DIR}}`: the entry file `{{ENTRY}}`, which exports the function below, any other
`.ts` files it imports from the same folder, and a `.css` file imported from the code (`import './styles.css';`)
if the feature needs styles of its own, as a plugin bundles them. Write nowhere else.

To learn the library, start with u2's `AGENTS.md`, then read its `README.md`, `docs/recipes/*.md`, `manifest.json`
and `src/dg/README.md`, all in `{{PUBLIC}}/libraries/u2/`; and the source of u2 (`{{PUBLIC}}/libraries/u2/src/`,
`{{PUBLIC}}/libraries/u2/css/`) and of the JS API (`{{PUBLIC}}/js-api/`). Read nothing else: no other plugin or
library, nothing under `libraries/u2/tests/` or `libraries/u2/eval/`, no git history.

The plugin's tests exercise the feature headlessly, in a page that stands in for the platform. The page has u2 and
these platform APIs, and no others (no other server calls, no viewers, no `ui.dialog` or `ui.input`):

- `DG.DataFrame` (`name`, `rowCount`, `columns.toList()`, `names()`, `byName(name)`, `bySemType(semType)`,
  `col(name)`, `rows.addNew(values, notify)`, its `on…` events), `DG.Column` (`name`, `type`, `semType`,
  `get(i)`), `DG.FileInfo` (`fullPath`, `fileName`, `DG.FileInfo.fromBytes`), `DG.Func` (`name`, `friendlyName`,
  `inputs`, `prepare(values)`), `DG.FuncCall` (`inputs`, `setParamValue`), `DG.Property` (`name`,
  `propertyType`, `defaultValue`, `caption`), `DG.TableView` (`dataFrame`, `getRibbonPanels`, `setRibbonPanels`,
  `close()`), `DG.View`, `DG.ViewBase`, `DG.TYPE`, `DG.SEMTYPE`;
- `grok.shell` (`v`, `views`, `tables`, `tableNames`, `table(name)`, `addTable`, `info`, `warning`, `error`),
  `grok.events` (`onCurrentViewChanged`, `onViewAdded`, `onViewRemoved`, `onTableAdded`, `onTableRemoved`),
  `grok.dapi.files` (`exists`, `list`, `readAsText`), `grok.chem.drawMolecule(smiles, width, height)`;
- `ui.pickFileFromFiles()`, the platform's Datagrok Files browser.

The platform hosts what the plugin's function returns or shows, and disposes it (`dispose()`) when it goes away.

Nothing here runs the plugin, so read the sources when you need to be sure of an API. Do not ask questions: where
the request leaves something open, decide as the platform's own UI would. When the code is written, reply with one
short paragraph on what you built.

{{TASK}}

---

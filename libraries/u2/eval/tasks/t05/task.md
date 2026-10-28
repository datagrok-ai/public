# Run a function with its settings

The compute pane runs a function on one of the tables the user has open. Show a dialog with the function's
settings, one field per input parameter, and **RUN**. The pane's function takes a `table` (a dataframe, one of
the open tables) and a `threshold` (a number, 0.5 by default). RUN calls `onRun(call)` with a `DG.FuncCall` of
`func` holding the user's choices: the table as the `DG.DataFrame` itself, and the threshold.

The plugin calls this, and it shows the dialog and returns it:

```ts
export function settingsDialog(func: DG.Func, onRun: (call: DG.FuncCall) => unknown): Dialog
```

# Compute settings

A campaign's compute settings say which functions run over the campaign's table, and with which parameters.
`store.load()` fetches the functions available and the saved settings, which takes a moment, so the pane shows that
it is loading meanwhile. Then it lists the functions:

- Each function is a check captioned with its friendly name, on for the functions in the saved settings, with its
  parameters under it: one field per input parameter, captioned with the parameter's name, except the first
  dataframe, which is always the campaign's `table`. Another dataframe parameter is one of the tables open in the
  session. A function's parameters can be changed only while it is checked.
- **SAVE** calls `store.save(settings)` with the checked functions' parameter values, by function and parameter
  name, a table by its name: `{descriptors: {threshold: 0.5}}`. It is off while no function is checked.
- The pane is at most 480px wide, and SAVE sits under the functions, at the right.

The plugin calls:

```ts
type Settings = {[func: string]: {[param: string]: unknown}};

interface ComputeStore {
  load(): Promise<{functions: DG.Func[], saved: Settings}>;
  save(settings: Settings): Promise<void>;
}

export function computeSettings(table: DG.DataFrame, store: ComputeStore): Control
```

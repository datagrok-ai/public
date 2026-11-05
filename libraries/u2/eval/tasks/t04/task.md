# Which column holds the structures

Before a merge, the user confirms which column of the campaign's table holds the structures. Return an input
captioned **Molecules** whose value is the column's name. It offers the table's molecule columns and starts on
the first one it offers, so the user usually just accepts it.

The plugin calls:

```ts
export function moleculeColumnInput(table: DG.DataFrame): Input<string | null>
```

---
title: "Register context actions"
---

To register a context-specific action that would get offered to a user when right-clicking
on the item or expanding the **Actions** pane, add the following tags to the function:

* `meta.action`: Text to be shown
* `input`: Has to be exactly one input, annotated with `semType`.

The following example registers a **Use as filter** action for molecules:

```js
//name: Use as filter
//description: Adds this structure as a substructure filter
//meta.action: Use as filter
//input: string mol { semType: Molecule }
export function useAsSubstructureFilter(mol: string): void {
  let tv = grok.shell.tv;
  let molCol = tv.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
  tv.getFiltersGroup({createDefaultFilters: false}).add({
    type: FILTER_TYPE.SUBSTRUCTURE,
    column: molCol.name,
    columnName: molCol.name,
    molBlock: molToMolblock(mol, getRdKitModule())
  });
}
```

This is the end result (note the **Use as Filter** action in the right panel):

![The custom context action for cells with molecules](context-actions.png)

See also:

* [Datagrok JavaScript development](../develop.md)
* [Test packages](test-packages.md)

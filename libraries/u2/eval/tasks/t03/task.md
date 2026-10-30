# Merge a file into the campaign

The **Merge with table** dialog asks for the file to merge into the campaign: one field, **Source file**, and
**MERGE**.

- Users keep these files on their computer or in Datagrok Files, and often paste a file share path such as
  `System:DemoFiles/demog.csv`. The field takes a file from any of the three.
- MERGE calls `onMerge(file)` with the chosen `DG.FileInfo`, and is off while no file is chosen. A file from
  Datagrok Files is passed with its share path (`fullPath`), so the merge can be repeated later.

The plugin calls this, and it shows the dialog and returns it:

```ts
export function mergeDialog(onMerge: (file: DG.FileInfo) => unknown): Dialog
```

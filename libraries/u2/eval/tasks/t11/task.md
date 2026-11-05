# Template editor

The editor of a new template, as the designer drew it:

- Three fields: **Name**, **Key** and **Description**. Name and Description take the whole width of their row; Key
  stays short, about 120px, since it holds a few capital letters.
- The editor is never wider than 640px.
- **CANCEL** and **CREATE** sit at the bottom, right-aligned, CANCEL first, 8px apart. CREATE calls
  `onCreate({name, key, description})`, CANCEL calls `onCancel()`.

The plugin calls:

```ts
export function templateEditor(onCreate: (values: {name: string, key: string, description: string}) => unknown,
  onCancel: () => void): HTMLElement | Control
```

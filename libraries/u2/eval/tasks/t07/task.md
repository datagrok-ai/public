# New template

A form for a new template: **Name** (required) and **Key** (required, 2 to 4 capital letters such as `KHL`).
**SAVE** calls `onSave({name, key})`; it stays off while the form cannot be saved.

The plugin calls:

```ts
export function newTemplateForm(onSave: (values: {name: string, key: string}) => unknown): HTMLElement | Control
```

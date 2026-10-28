# Campaigns on the app's home

The app's home lists the user's campaigns: a table with the columns Name, Code, Status and Created, grouped by
Status, each group headed by its status and how many campaigns it holds. As in the platform's other lists, the row
under the pointer is highlighted.

The plugin calls:

```ts
export function campaignsTable(campaigns: {code: string, name: string, status: string, created: Date}[]):
  HTMLElement | Control
```

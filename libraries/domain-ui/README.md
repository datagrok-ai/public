# @datagrok-libraries/domain-ui

Composable UI over [entity-mapped domain tables](https://datagrok.ai/help/datagrok/) —
the `grok.dapi.domains` data plane's UI counterpart.

Everything is reflective: components take columns, labels, choices, validation rules and
permissions from the runtime domain registry, so they work on any registered table with no
codegen and no hand-wired permission checks.

## Install

```
npm install @datagrok-libraries/domain-ui
```

## An editable grid, saving as one transaction

```ts
import {DomainGrid} from '@datagrok-libraries/domain-ui';

const grid = await DomainGrid.create(grok.dapi.domains.table('grit.issue'));
grok.shell.newView('Issues', [grid.root]);
```

Edited cells highlight amber, invalid ones red (with the server's own message in the
tooltip), rows marked for deletion disappear from the view but stay undoable, and Save
writes the whole batch as ONE `/transaction` — audit rows share a `tx_id`, a version
conflict goes through the platform's standard reload/overwrite dialog.

## The state model

`DomainFrameEditor` is the single writer of three service columns it attaches to the frame:

| Column | Content |
|---|---|
| `~state` | `'' \| 'new' \| 'modified' \| 'deleted'` per row |
| `~changes` | JSON — the ORIGINAL value of every changed cell (sparse) |
| `~errors` | JSON — per-cell validation / conflict messages |

Each is tagged out of binary AND CSV export, so the editing state is memory-only: a saved
project, `toByteArray()`, `toCsv()`, an export or a `batch()` upload built from the frame
never carry it.

Write through `setValue()` (programmatic) or `beginEdit()` + `commitEdit()` (what a grid
does). Writing a cell directly on the DataFrame bypasses tracking and is silently not saved.

## Refreshing discards pending edits — by design

`refresh()` re-runs the query and rebuilds the frame and its state from scratch. There is no
merge. Deciding whether it is safe to refresh is the caller's job:

```ts
if (editor.isDirty && !await confirmDiscard())
  return;
await editor.refresh();
```

## Permissions

Every affordance comes from `DG.DomainTableCapabilities`: no `canEdit` gives a read-only
grid, no `canInsert` removes Add row, no `canDelete` removes Delete row, and columns outside
`writableColumns` stay read-only and never appear in a payload.

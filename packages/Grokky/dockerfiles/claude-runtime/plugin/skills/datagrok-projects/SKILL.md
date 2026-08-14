---
name: datagrok-projects
description: Save, publish, create, attach to, or share a Datagrok project, including turning data sync on. The save and share flows drive the platform dialogs through their AI functions.
---

# datagrok-projects

## Save / publish — drive the platform's Save dialog

The Save Project dialog is the platform's save/publish entry and the only reliable way to get
per-table data sync, creation-script dependency linking, and the correct upload order. Do not
hand-write project/table saves with `grok.dapi` — open the dialog and drive it through its
functions.

1. **Open the dialog** (`datagrok_exec`) — do NOT await the dialog promise (it resolves only when
   the dialog closes):

   ```js
   const project = grok.shell.project;
   DG.Project.showSaveDialog({
     tables: [t], views: [view],                              // add more pairs for multi-table saves
     project: project.isOnServer ? project : undefined,       // re-publish into the SAME project
     name: project.isOnServer ? undefined : '<name the user gave>',
   });
   return {dialogOpened: true};
   ```

2. **Find it**: `list_view_widgets` — the dialog is an extra root ref `dlg0` ("Save project");
   note also the child ref of its entity-list widget (type `ProjectEntityMoveWidget`).

3. **Read its state**: `get_view_function_result(name: 'getProjectSaveInfo', widget: 'dlg0')` —
   name (with validity), description, presentation mode, and every entity to be saved, with
   `dataSync` / `dataSyncPossible` per table.

4. **Adjust only what the user asked for**, via `call_view_function`:
   - on `dlg0`: `setProjectName`, `setProjectDescription`, `setPresentationMode`, `setSaveMode`
     (when the dialog offers it)
   - on the entity-list widget ref: `listEntities`, `setDataSync(table, enabled)`,
     `setEntityAction(name, action)`

   Every function returns `{success: false, error: ...}` with the reason when something can't be
   done (e.g. data sync without a creation script) — report that reason; never work around it
   with raw tag/meta writes.

5. **Confirm**: `call_view_function(name: 'clickButton', parameters: {name: 'OK'}, widget: 'dlg0')`
   — applies the save and closes the dialog. If the user wanted to review first, stop before this
   step and say the dialog is ready.

6. **Verify server-side** (`datagrok_verify`) — the upload is asynchronous, so if the first
   re-read misses, confirm the dialog has closed and retry once:

   ```js
   const p = await grok.dapi.projects.filter('friendlyName = "<name>"').first();
   return p != null;
   ```

   When data sync was requested, also verify the table:
   `(await grok.dapi.tables.find('<tableInfoId>')).tags['.data-sync'] === 'sync'`.

## Data sync — the one precondition

A synced table re-runs its **creation script** when the project opens. The platform records that
script only when the query/script runs through its own pipeline — `grok.functions.call` /
`grok.data.query` run as processed calls and record nothing. So for "query the database and save
it with data sync":

1. Make sure the query is **saved on the server**: create it with
   `connection.query('<name>', '<sql>')`, `q = await grok.dapi.queries.save(q)`, then re-read by
   **id** (`grok.dapi.queries.find(q.id)`) — never by name, duplicates shadow name lookups.
2. Run the server copy so the result opens in the workspace **with** a creation script:

   ```js
   const call = q.prepare({});
   await call.call(false, undefined, {processed: false});
   ```

3. Open the Save dialog (above) and `setDataSync(table, true)`. If it answers
   "no creation script", step 2 didn't happen — fix that instead of forcing tags.

## Create / attach / share

- **Create** an empty project: `datagrok_projects(op: "create", args: {name})` once.
- **Attach** the current table/view to a project: use the Save dialog flow above.
- **Share**: `datagrok_projects(op: "share", args: {projectId, groups, access})` —
  `access` `"View"`|`"Edit"`, default `"View"`; ask when no group is named.
  When the platform's **Share dialog** is open instead (a `dlg<N>` ref), drive it the same way:
  `getShareInfo`, `addShareGrantee(name)`, `setShareAccess(access)`, `setShareMessage(message)`,
  then `clickButton OK`.

## Don't

- Hand-write saves with `grok.dapi.tables.save` / `uploadDataFrame` / `projects.save` — the
  dialog handles upload order, table meta, and sync dependencies; manual writes corrupt them.
- Set `.data-sync` / `.script` tags directly — use `setDataSync` and let it refuse when the
  precondition fails.
- Re-run a save after a timeout — re-read state first (`getProjectSaveInfo`, or check the dialog
  closed); a second dialog or a duplicate project is worse than a late answer.
- Promise data sync before `setDataSync` (or `dataSyncPossible`) confirms it.

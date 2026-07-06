---
name: datagrok-projects
description: Save, create, attach to, or share a Datagrok project. Each step has a trigger phrase; run only the triggered step, then stop.
---

# datagrok-projects

Run only the triggered step. Never chain ahead or infer the next step.

## Routing

| Trigger                                                | Run             |
|--------------------------------------------------------|-----------------|
| "save the/my project / save my work" (no "new")        | Step 1          |
| explicit "create / make / new project"                 | Step 2          |
| "attach the (current) table / view"                    | Step 3          |
| "share the project / with <group>"                     | Step 4          |
| Step 1 returned `onServer: false`                      | Step 2 + Step 3 |

"Save" never means create — a new project on a "save" request is a clone; the original stays unsaved.

Allowed chains only: Step 1 → Step 2 + Step 3 on `onServer: false`; Step 2 → Step 3 when the
same message says "with the current table/view/dataframe". Otherwise stop after the step.

## Step 1 — Save (update the open project)

`datagrok_exec`:

```js
// Save the LIVE workspace project — a fresh DG.Project.create() would be a clone
// (public repo: UsageAnalysis TestTrack/Charts/radar-save-reopen-bug-spec.ts)
const project = grok.shell.project;
if (!project.isOnServer)
  return {onServer: false};
// Persistence sequence per Bio/src/tests/projects-tests.ts saveAndOpenProject()
const tableInfo = t.getTableInfo();
const viewInfo = view.getInfo();
await grok.dapi.tables.uploadDataFrame(t);
await grok.dapi.tables.save(tableInfo);
await grok.dapi.views.save(viewInfo);
await grok.dapi.projects.save(project);
return {projectId: project.id, name: project.friendlyName, updated: true};
```

- The entities keep their server ids → each save updates the same project (UI SAVE button behavior).
- Never `DG.Project.create()` or `create_project` here; no `addChild` — the entities already belong to the project.

Result:
- `updated: true` → show the project with `datagrok_show_entities`, confirm saved. Stop.
- `onServer: false` → workspace is the local scratchpad; nothing to update. Ask for a project
  name if the user gave none, then run Step 2 (create) + Step 3 (attach).

## Step 2 — Create

`create_project(name)` once → keep returned id as `projectId`. Stop (unless chaining per Routing).

## Step 3 — Attach current table + view

Requires `projectId`. `datagrok_exec` (substitute `projectId`):

```js
const project = await grok.dapi.projects.find('<projectId>');
const tableInfo = t.getTableInfo();
await grok.dapi.tables.uploadDataFrame(t);
await grok.dapi.tables.save(tableInfo);
const viewInfo = view.getInfo();
await grok.dapi.views.save(viewInfo);
project.addChild(tableInfo);
project.addChild(viewInfo);
await grok.dapi.projects.save(project);
return {projectId: project.id, ok: true};
```

## Step 4 — Share

`share_project(projectId, groups, access)` — `access` `"View"`|`"Edit"`, default `"View"`. No group named → ask.

## Don't

- Run `create_project` on "save" while a server project is open — update in place (Step 1).
- `view.saveLayout()` / `grok.dapi.layouts.save()` — reopens as bare grid.
- Create or share in `datagrok_exec` — use MCP tools.
- Pre-announce later steps; call `create_project` twice.
- On MCP tool error: surface verbatim, stop. No `datagrok_exec` fallback.

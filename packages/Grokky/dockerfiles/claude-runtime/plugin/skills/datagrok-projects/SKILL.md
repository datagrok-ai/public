---
name: datagrok-projects
description: Create, attach to, or share a Datagrok project. Each step has a trigger phrase; run only triggered steps, in order, then stop.
---

# datagrok-projects

Run only the triggered step. Never chain ahead or infer the next step.

## Step 1 — Create

Trigger: "create / make / new project".

`create_project(name)` once → keep returned id as `projectId`. Stop.

Chain into Step 2 only if the same message says "with the current table/view/dataframe".

## Step 2 — Attach current table + view

Trigger: "attach the (current) table / view"; `projectId` known.

`datagrok_exec` (substitute `projectId`):

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

All four saves required. Omit `uploadDataFrame` → "insufficient privileges". Omit `views.save` → reopens as bare grid (no viewers). Table only → drop the two `viewInfo` lines.

## Step 3 — Share

Trigger: "share the project / with <group>".

`share_project(projectId, groups, access)` — `access` `"View"`|`"Edit"`, default `"View"`. No group named → ask.

## Don't

- `view.saveLayout()` / `grok.dapi.layouts.save()` → bare grid.
- Create or share in `datagrok_exec` — use MCP tools.
- Pre-announce later steps; call `create_project` twice.
- On MCP tool error: surface verbatim, stop. No `datagrok_exec` fallback.

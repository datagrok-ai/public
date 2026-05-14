---
name: datagrok-projects
description: Use whenever the user asks to create, assemble, populate, or share a Datagrok project — especially "create a project with the current table / view / layout". Walks through the exact, serialized recipe so MCP tools and the single datagrok-exec block fire in the right order.
---

# datagrok-projects

## Recipe (strict order — no parallelization)

### Step 1 — create the project

Call `create_project(name)` EXACTLY ONCE. Save the returned `id` as `projectId`.

Do NOT call `create_project` twice. Do NOT fire it in parallel with anything else.
Do NOT emit a datagrok-exec block in the same turn as `create_project` (you do not
yet have `projectId`, and the exec block's return value is not visible until the
next turn anyway).

### Step 2 — (only if the user said "with the current table / view / layout / dataframe")

Emit ONE `datagrok-exec` block that saves the layout and returns the ids:

```datagrok-exec
const layout = view.saveLayout();
await grok.dapi.layouts.save(layout);
return {tableInfoId: t.getTableInfo().id, layoutId: layout.id};
```

The block's return value comes back to you on the NEXT turn as a fenced
`<exec-result>` JSON object. Read `tableInfoId` and `layoutId` from there.

You do NOT need to ask the user for the ids — they arrive automatically.

### Step 3 — attach the entities

Call `add_entity_to_project(projectId, tableInfoId)` AND
`add_entity_to_project(projectId, layoutId)`. Two calls — one per entity.

### Step 4 — share (if the user asked)

Call `share_entity(projectId, groupName, access)` where:
- `projectId` is the UUID from step 1 (preferred over `"Namespace:Name"` form),
- `groupName` is the plain group name the user mentioned (e.g. `"toxicology-review"`),
- `access` is `"View"` or `"Edit"`.

Do NOT call `get_group` first — `share_entity` accepts plain names. Only call
`get_group` if `share_entity` throws a group-not-found error and you need to
verify the spelling.

### Step 5 — return the link

Call `get_project_link(projectId)` ONCE. Include the returned `url` in your reply.

## Forbidden patterns

- `grok.dapi.projects.save(project)` inside datagrok-exec — bypasses governance.
- `project.addChild(...)` / `project.addLink(...)` inside datagrok-exec — bypasses MCP.
- `grok.dapi.permissions.grant(...)` inside datagrok-exec — bypasses share governance.
- Calling `create_project` more than once in the same workflow.
- Firing `create_project` in parallel with the datagrok-exec block or with
  `add_entity_to_project`. Each step must finish, and its result must be visible,
  before the next begins.
- Calling `get_group` defensively before `share_entity`.
- Asking the user for the table or layout id — the datagrok-exec block returns them.

## Failure handling

If any MCP tool throws, surface the error to the user verbatim and stop. Do NOT
silently fall back to a datagrok-exec block that does the same operation via the
in-browser API.

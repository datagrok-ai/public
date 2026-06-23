---
name: datagrok-projects
description: Use whenever the user asks to create, attach to, populate, or share a Datagrok project. Each step below has its OWN trigger phrase — do only the steps the user actually asked for, in order, then stop and wait for the next instruction. Never chain ahead.
---

# datagrok-projects

## Golden rule

Each step is gated by a SPECIFIC trigger phrase. Run only the steps the user
explicitly asked for. **Do not infer next steps from recent conversation
context** (e.g. earlier mentions of a table, a previous attach run, a layout
that's open). When in doubt, do less and let the user ask for the next step.

## Step 1 — Create

**Trigger:** the user asks to "create a project" / "make a project" / "new project".

**Do:** call `create_project(name)` ONCE. Save the returned `id` as `projectId`.

**Stop here.** Do NOT pre-announce attach or share. Do NOT emit a
`datagrok-exec` block. Do NOT call `attach_entity_to_project`. Wait for the
user's next message.

**Exception:** Step 1 may chain DIRECTLY into Step 2 ONLY when the user's same
message explicitly includes "with the current table", "with the current view",
"with the current layout", "with this dataframe", or an equivalent attach
phrase. A bare "create a project Foo" is Step 1 only.

## Step 2 — Persist live entities to the server

**Trigger:** the user asks to "attach the (current) table / layout / view" OR
the create message included one of the explicit attach phrases above. AND
`projectId` from Step 1 is known.

**Do:** call `datagrok_exec` with code that persists the live entities and
returns their server-side ids. `t.getTableInfo().id` is a client-side stub —
attaching it will fail with "insufficient privileges to edit". You MUST call
`uploadDataFrame` first.

```js
const tableInfoId = await grok.dapi.tables.uploadDataFrame(t);
const layout = view.saveLayout();
await grok.dapi.layouts.save(layout);
return {tableInfoId, layoutId: layout.id};
```

The ids come back synchronously in the tool result as `returnValue`.
Do NOT ask the user for the ids.

If the user only asked about the table (not the layout), drop the layout lines.
If only about the layout, drop the table lines.

## Step 3 — Attach

**Trigger:** Step 2 just completed and you have its returned ids.

**Do:** for each id you persisted in Step 2, call
`attach_entity_to_project(projectId, <id>)`. One call per entity.

Stop after the attach calls. Do NOT volunteer to share.

## Step 4 — Share

**Trigger:** the user explicitly asks to "share the project" / "share with
<group>" / "give <group> access".

**Do:** call `share_project(projectId, groups, access)` where:
- `projectId` is the UUID from Step 1,
- `groups` is the array of plain group names the user named,
- `access` is `"View"` or `"Edit"` (defaults to `"View"`).

If the user said "share the project" without naming a group, ask them which
group(s) before calling.

## Forbidden patterns

- Pre-announcing future steps in the same turn ("I'll create the project, then
  attach the table, then share..."). Announce only the step you're doing now.
- Chaining Step 2/3/4 from a bare "create a project Foo" message.
- Returning `t.getTableInfo().id` from the datagrok_exec code — client stub.
  Always use `await grok.dapi.tables.uploadDataFrame(t)`.
- `grok.dapi.projects.save(project)` in datagrok_exec code — bypasses governance.
- `project.addChild(...)` / `project.addLink(...)` in datagrok_exec code — bypasses MCP.
- `grok.dapi.permissions.grant(...)` in datagrok_exec code — bypasses share governance.
- Calling `create_project` more than once in the same workflow.
- Asking the user for the table or layout id — the tool result contains them.

## Failure handling

If any MCP tool throws, surface the error to the user verbatim and stop. Do NOT
silently fall back to a datagrok-exec block that does the same operation via the
in-browser API.

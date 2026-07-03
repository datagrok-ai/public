---
name: datagrok-entities
description: Call datagrok_show_entities after any MCP tool result that returns entity data (files, scripts, queries, connections, projects, spaces, groups, users). Always use the tool — never list entities as plain text, markdown links, or bullet lists.
---

# datagrok-entities

When an MCP tool result contains entity data, call `datagrok_show_entities` immediately.
The tool renders interactive cards the user can click to open the entity — always prefer
it over prose.

## Tool signature

```
datagrok_show_entities(entities: EntityRef[])
```

## Supported types

| type         | required                    | optional                      |
|--------------|-----------------------------|-------------------------------|
| `file`       | `connector`, `path`, `name` | `isDirectory`                 |
| `script`     | `id`, `name`                | —                             |
| `query`      | `id`, `name`                | —                             |
| `connection` | `id`, `name`                | —                             |
| `project`    | `id`, `name`                | —                             |
| `space`      | `id`, `name`                | —                             |
| `group`      | `id`, `name`                | —                             |
| `user`       | `id`, `name`                | —                             |

## Notes

- **`#type` passthrough**: MCP results return entities with a `#type` field in PascalCase (e.g. `"FileInfo"`, `"DataQuery"`, `"UserGroup"`). Pass it as-is in the `#type` field — the tool normalizes it to the correct `type` value automatically. You may also pass the normalized `type` string directly if you already know it.
- Pass `id` and all other identifying values **exactly** from the MCP response — never invent.
- For files, `connector` is the connection name string (e.g. `"System:DemoFiles"`) and `path` is the file path relative to the connector root.
- For users, use `friendlyName` as the `name` field — users have no `name` property in the MCP response, only `friendlyName` and `login`.
- A single call may mix entity types — the renderer handles them.
- After calling the tool, you may add a short prose summary if useful, but the tool call must come first.

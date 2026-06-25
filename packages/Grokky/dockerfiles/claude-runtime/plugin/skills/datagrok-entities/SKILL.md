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
| `file`       | `connector`, `path`, `name` | `isDirectory`, `size`         |
| `script`     | `id`, `name`                | `language`                    |
| `query`      | `id`, `name`                | `connectionName`              |
| `connection` | `id`, `name`                | `dataSource`                  |
| `project`    | `id`, `name`                | —                             |
| `space`      | `id`, `name`                | —                             |
| `group`      | `id`, `name`                | —                             |
| `user`       | `id`, `name`                | —                             |

## Notes

- Pass `id` and all other identifying values **exactly** from the MCP response — never invent.
- For files, `connector` is the connection name (e.g. `System:DemoFiles`) and `path` is
  relative to the connector root.
- A single call may mix entity types — the renderer handles them.
- After calling the tool, you may add a short prose summary if useful, but the tool call
  must come first.

---
name: datagrok-entities
description: Use after a Datagrok MCP tool result contains entity data — files, scripts, queries, connections, projects, or spaces. Wrap those entities in a datagrok-entities fenced block so they render as interactive cards instead of markdown links or bullet lists. Required whenever the Datagrok MCP returns entities; never paraphrase them as plain text.
---

# datagrok-entities

Emit a `datagrok-entities` fenced block with a JSON array. Each entry renders
as an interactive card the user can click to open the entity.

## Block format

```datagrok-entities
[{"type":"file","connector":"System:DemoFiles","path":"datasets/demog.csv","name":"demog.csv","size":12345}]
```

## Supported types

| type         | required               | optional             |
|--------------|------------------------|----------------------|
| `file`       | `connector`, `path`, `name` | `isDirectory`, `size` |
| `script`     | `id`, `name`           | `language`           |
| `query`      | `id`, `name`           | `connectionName`     |
| `connection` | `id`, `name`           | `dataSource`         |
| `project`    | `id`, `name`           | —                    |
| `space`      | `id`, `name`           | —                    |

## Notes

- `id` and other identifying values come from the MCP tool response — pass them
  through exactly, never invent.
- For files, `connector` is the connection name (e.g. `System:DemoFiles`) and
  `path` is relative to the connector.
- A single block may contain mixed entity types — the renderer handles them.

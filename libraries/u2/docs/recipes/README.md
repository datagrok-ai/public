# u2 recipes

Canonical shapes for building Datagrok UI with u2 — for humans and for AI assistants. Each
recipe is one file, so it can be referenced by name. Read the recipe before composing the
corresponding UI; the components' own reference lives in `manifest.json` (each entry's
`usage` field carries per-component judgment, these files carry the composition-level kind).

The governing principle: **u2 apps run inside the real Datagrok shell. Don't recreate the
shell's chrome inside your view — use it.**

| Recipe | When |
|---|---|
| [shell-integration](shell-integration.md) | Any app or view meant to run in the platform |
| [entity-browser](entity-browser.md) | A list of server objects with details — the master-detail shape |
| [list-item-rendering](list-item-rendering.md) | Rendering rich items (reports, mail, issues) in a list |
| [catalog-table](catalog-table.md) | A catalog / records table: campaigns or profiles to sort, group, choose columns of and act on |
| [forms-and-details](forms-and-details.md) | Editing an object, or showing a read-only detail pane; multi-step flows |
| [specs](specs.md) | Writing a dg-ui/1 spec: live binds, structural props, events, plugins (`installEms`) |

Whole, copy-ready code for the everyday jobs: [docs/samples/](../samples/README.md). Every TypeScript block here
and every sample is compiled against u2's sources and the live ones run under the testing kit
(`tests/recipes.test.js`), so what these files show is what u2 does.

Domain (EMS) recipes: [libraries/ems/docs/recipes/](../../../ems/docs/recipes/).

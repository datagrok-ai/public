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
| [forms-and-details](forms-and-details.md) | Editing an object, or showing a read-only detail pane |
| [crud-app](crud-app.md) | A domain (EMS) table as an app with no code — `schema.json` + `table.app()`: what the schema declares, what the app gives |
| [spec-app](spec-app.md) | The same app as a designer-editable `dg-ui/1` spec — `u2-domain-source` + the `u2-domain-*` tags, one session per page |
| [hierarchies](hierarchies.md) | A self-referencing table as a tree — `hierarchy: true`, `pathTo`, the `under` subtree filter, `domains.tree` driving a collection |
| [custom-app](custom-app.md) | The app in code — actions, validators, renderers, presets, shortcuts, a `DomainApp` subclass, master–detail in one session |

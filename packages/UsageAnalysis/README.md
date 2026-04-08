# Usage analysis

Usage analysis is a package that shows statistics of Datagrok platform usage.

In the different views, you can find different information about the platform events and errors.

## Clicks

The **Clicks** tab aggregates user interaction events (clicks, inputs, commands, dialogs) broken down by their full UI path (e.g. `Inspector / Click Events / Date range`). Useful for identifying heavily-used and under-used UI areas.

The path is decomposed into **Level 1–5** columns for hierarchical filtering, and a **Source** column distinguishes platform-built-in actions from package-contributed ones.

Filters on the left:
- **Min count** slider — hides low-frequency rows (noise reduction)
- **event_type** — filter by click, input, command, dialog show/close/ok, etc.
- **Source** — Platform vs Package
- **~Package Name** — filter by specific package
- **Level 1–3** — drill into specific UI areas

The tree map at the bottom shows click frequency by Level 1–3, with cell size proportional to count.

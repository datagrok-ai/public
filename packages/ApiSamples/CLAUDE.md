# ApiSamples

ApiSamples is the canonical reference package for Datagrok API usage. It contains ~374 standalone JavaScript example
scripts organized by API area, plus a test harness that runs each script automatically against a live Datagrok instance.

Direct link to open samples tree: `${GROK_URL}/js?browse=samples`.

## Adding a New Script

Scripts are plain `.js` files dropped into the appropriate subfolder under `scripts/`. They execute as top-level async
code when invoked by the platform. 

Keep them very short and clean.

Optional first line `//api: DG.Viewer.fromType, ui.input.string` names the members the sample demonstrates
(qualified names, comma-separated). The JS API inventory (`core/docs/reviews/js-api-audit/data/inventory.cjs`)
builds its samples-per-member index from that line, and falls back to the qualified `DG.`/`ui.`/`grok.` usages in
the source when it is absent.

### Script categories

| Folder         | What belongs here                                                   |
|----------------|---------------------------------------------------------------------|
| `dapi/`        | Server API (users, projects, files, layouts, messaging)             |
| `dapi/domains/`| Domain schemas: row CRUD, queries, the domain UI and app framework  |
| `data-frame/`  | DataFrame construction, modification, filtering, aggregation, joins |
| `data-access/` | Database and external API access                                    |
| `ui/`          | UI components (buttons, inputs, dialogs, viewers, layouts)          |
| `grid/`        | Table grid customization                                            |
| `events/`      | Event subscriptions                                                 |
| `functions/`   | Function API                                                        |
| `domains/`     | Bio, chem, data-science domain APIs                                 |
| `scripting/`   | Script parameter/input/output patterns                              |
| `performance/` | Benchmarks (large data)                                             |
| `shell/`       | Shell view management, notifications                                |
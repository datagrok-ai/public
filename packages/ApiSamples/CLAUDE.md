# ApiSamples

ApiSamples is the canonical reference package for Datagrok API usage. It contains ~374 standalone JavaScript example
scripts organized by API area, plus a test harness that runs each script automatically against a live Datagrok instance.

Direct link to open samples tree: `${GROK_URL}/js?browse=samples`.

## Adding a New Script

Scripts are plain `.js` files dropped into the appropriate subfolder under `scripts/`. They execute as top-level async
code when invoked by the platform. 

Keep them very short and clean.

### Script categories

| Folder         | What belongs here                                                   |
|----------------|---------------------------------------------------------------------|
| `dapi/`        | Server API (users, projects, files, layouts, messaging)             |
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
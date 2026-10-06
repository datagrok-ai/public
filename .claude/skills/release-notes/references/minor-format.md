# Minor release format (1.Y.0, for example 1.28.0)

## Classification: which section a ticket belongs in

Determine section by priority:

1. **Commit message prefix** — match against known prefixes:
   - `Grid:` → Grid
   - `Scatterplot:` / `Scatter plot:` → Scatterplot
   - `Filter:` / `Filters:` → Filter Panel
   - `Browse:` → Browse
   - `JS API:` → JS API
   - `Data Access:` / `Visual Query:` / `Connections:` / `Queries:` → Data Access
   - Package names (`Bio:`, `Chem:`, `Charts:`, `DiffStudio:`, `EDA:`, `PowerGrid:`,
     `UsageAnalysis:`, etc.) → corresponding Package section
   - Viewer-specific (`Box plot:`, `Line Chart:`, `PC plot:`, `Trellis:`, `Pivot table:`,
     `Tile Viewer:`, `Calendar:`, etc.) → viewer-specific or general Viewers section
   - Everything else platform-related → Platform

2. **Ticket summary prefix** — same prefix-matching logic applied to the issue title.

3. **Commit file paths** — run `git show --stat <hash>` for ambiguous commits:
   - `public/packages/<Name>/` → Package `<Name>`
   - `public/js-api/` → JS API
   - `core/client/d4/viewers/grid/` → Grid
   - `core/client/d4/viewers/scatter_plot/` → Scatterplot
   - `core/client/d4/viewers/filters/` → Filter Panel
   - `core/client/xamgle/web/browse/` → Browse
   - Other `core/client/` or `core/server/` → Platform

**Bug vs Feature classification:**
- Jira: `issuetype.name` == "Bug" → `* Fixed:` subsection
- GitHub: title starts with "Fix" / contains "Fixed" / has label "bug"
- Commit: message contains "Fixed" / "Fix" prefix

## Main updates section (Jira label `main`)

- Tickets with label `main` go into the "Main updates" section AND their respective category section.
  If none is labeled, propose the largest user-facing items of the draft.
- **Ordering**: group related items thematically (AI-related items together,
  then infrastructure, then UX). Most impactful / user-visible items first.
  Present the proposed order to the user for confirmation.
- **Selection criteria**: only strategic, cross-platform, user-facing features belong here.
  Viewer-specific improvements (even if labeled `main`) go in their viewer section only —
  do not promote them to Main updates unless they represent a platform-wide capability.

**Writing style for Main updates** (overrides general reformulation rules):
- Write in **user-centric language**: "lets you", "gives you", "give administrators", etc.
  Never use passive construction ("Added", "Implemented", "Introduced") as the opening.
- Use a **clean feature name** as the bold title — no component prefix like
  "Grid: column-based coloring". Write just the concept: "Column-based coloring".
  Exception: keep the component if essential to identify the feature.
- After the bold title, **explain the benefit** — what the user can now do and why it matters.
- Use **em dash (`—`)** to append secondary context or benefit:
  `**Feature name** lets you do X — which improves Y`.
- If a related doc page exists, add `For details, see [link]` at the end of the item.
- **Do not include viewer-specific improvements** even if those tickets have label `main`.

Example:
```markdown
* **Roles** give administrators a structured way to define access levels and assign permissions
  to groups of users at once, simplifying access management across the platform
* **Script layouts** let you save and restore viewer configurations, column coloring, and styles
  alongside your script, so your analytical setup is always ready to reuse
* **Click tracking** gives platform administrators usage analytics based on user interactions across the UI
```

## Reformulation rules (general)

Source of text: the **ticket summary** (GitHub issue title or Jira summary), not the
commit message. Commit messages are only for classification and ticket extraction.

- Convert the ticket summary to a past-tense completed action.
- **Vary the leading verb** — never repeat "Added" on every line:
  - `Introduced` — brand-new features
  - `Improved` — enhancements to existing functionality
  - `Enabled` — unlocking a capability
  - `Exposed` — newly available API surface
  - `Extended X with Y` — expanding existing functionality
  - `Implemented` — technical constructs (rarely needed)
- **Avoid "Added the ability to"** — use `Enabled`, `Introduced`, or `Extended` instead.
- GitHub issues: `[#NNN](https://github.com/datagrok-ai/public/issues/NNN): ` + reformulated title.
- GitHub bug issues under `* Fixed:`: `[#NNN](link): ` + rewrite as positive outcome.
- Jira-only items: reformulated summary, no ticket reference.
- Commits without any ticket: reformulate the commit message as completed action.

**Bug descriptions under `* Fixed:` must NOT start with "Fixed":**
- Wrong: `* Fixed: \n  * Fixed correct state application`
- Right: `* Fixed: \n  * Correct state application`

**Bug fix phrasing — describe the positive outcome, not the problem:**
- Wrong: `Query renaming issues in Browse`
- Right: `Renaming the Query after it's been used in a Project now updates the script`

## Thematic sub-groups (Platform and JS API only)

Group items under nested bullet sub-headers instead of a flat list:

- Platform: `AI`, `Access control`, `Search`, `Scripting`, `Data`, `Infrastructure`
- JS API: `Views and UI`, `Inputs`, `Data`, `Other`

```markdown
### Platform

* Access control
  * Improved role permissions management with a new permissions editor
  * Introduced granular entity permissions
* Scripting
  * Introduced multi-selection choices for `list<string>` in all scripting languages
```

## Viewer structure rule

- The `### Viewers` section contains **cross-viewer improvements**.
- Each major viewer always gets its own `####` subsection under `### Viewers`.
- Minor viewers with ≤ 3 trivial items can be folded into `### Viewers` as inline mentions.
- Viewer sections use `####`, not `###`.

## Section pattern

```
### Section Name

* Feature/improvement item 1
* Feature/improvement item 2
* Fixed:
  * Bug fix item 1
  * Bug fix item 2
```

## Full document structure

```markdown
## <release-date> Datagrok <version> release

<1-2 sentence high-level summary. Keep it generic — do NOT list specific features.
  Example: "...introduces AI-powered features, brings improvements to visualization,
  strengthens access control and analytics, and delivers a range of usability
  enhancements across the platform and packages."
  Mention broad themes, not individual items.>

### Breaking changes

<Items with BREAKING keyword or detected breaking logic. Omit section if none.>

### Main updates

* **Bold title** lets you / gives you ... — benefit or context.

### Platform

* items...
* Fixed:
  * items...

#### Data Access

* items...
* Fixed:
  * items...

#### Browse

* items...
* Fixed:
  * items...

### [JS API](https://datagrok.ai/help/develop/packages/js-api)

* items...
* Fixed:
  * items...

### Viewers

* cross-viewer improvements (annotation regions, zoom state, column selectors, etc.)
* Fixed:
  * items...

#### [Grid](../../visualize/viewers/grid.md)
#### [Scatterplot](../../visualize/viewers/scatter-plot.md)
#### [Line Chart](../../visualize/viewers/line-chart.md)
#### [Bar chart](../../visualize/viewers/bar-chart.md)
#### [Histogram](../../visualize/viewers/histogram.md)
#### [Box plot](../../visualize/viewers/box-plot.md)
#### [Trellis plot](../../visualize/viewers/trellis-plot.md)
#### [Filter Panel](../../visualize/viewers/filters.md)

### Packages

#### [PackageName](https://github.com/datagrok-ai/public/tree/master/packages/PackageName)

* items compiled from submodule commits
```

## Community forum links

When a notable new viewer capability has a community.datagrok.ai post, add:
`see [FeatureName updates](https://community.datagrok.ai/t/...)` at the end of the item.
These links are provided by the user or found in commit messages — do not fabricate.

## Viewer help link reference

| Viewer | Path |
|---|---|
| Grid | `../../visualize/viewers/grid.md` |
| Scatterplot | `../../visualize/viewers/scatter-plot.md` |
| Filter Panel | `../../visualize/viewers/filters.md` |
| Box plot | `../../visualize/viewers/box-plot.md` |
| Line Chart | `../../visualize/viewers/line-chart.md` |
| PC plot | `../../visualize/viewers/pc-plot.md` |
| Trellis plot | `../../visualize/viewers/trellis-plot.md` |
| Pivot table | `../../visualize/viewers/pivot-table.md` |
| Tile Viewer | `../../visualize/viewers/tile-viewer.md` |
| Pie chart | `../../visualize/viewers/pie-chart.md` |

Package link pattern: `https://github.com/datagrok-ai/public/tree/master/packages/<Name>`

## Latest version Docker image table

The `## Latest version` table at the top of `release-history.md` is updated **at the very end**
of the release process, after Docker images are published. Format:

```markdown
| Service                                                                   | Docker Image                                                                                      |
|---------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------|
| [Datagrok](../../develop/under-the-hood/infrastructure.md#1-core-components) | [datagrok/datagrok:<version>](https://hub.docker.com/r/datagrok/datagrok)                         |
| [Grok Connect](../../develop/under-the-hood/infrastructure.md#3-external-database-connectivity) | [datagrok/grok_connect:<version>](https://hub.docker.com/r/datagrok/grok_connect)                 |
| Grok Spawner                                                              | [datagrok/grok_spawner:<version>](https://hub.docker.com/r/datagrok/grok_spawner)                 |
```

Fetch current versions when ready:
```bash
curl -s "https://hub.docker.com/v2/repositories/datagrok/datagrok/tags/?page_size=5&ordering=last_updated" | jq -r '.results[].name' | head -5
```

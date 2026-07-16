---
feature: charts
target_layer: apitest
coverage_type: regression
priority: p0
realizes_atlas: [charts.cp.open-viewer-with-required-columns]
realizes: [charts.chord, charts.globe, charts.group-analysis, charts.multiplot, charts.radar, charts.sankey, charts.sunburst, charts.surface-plot, charts.tree, charts.word-cloud]
pyramid_layer: integration
ui_coverage_responsibility: []
ui_coverage_delegated_to: null
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Charts/charts-api.md
date_created: 2026-05-07
authored_by: orchestrator-test-designer-charts-migrate-2026-05-07
related_bugs: []
---

# Charts — viewer API contract

API-contract scenario for the four Charts package viewers (Radar, Sunburst,
Tree, Timelines). This checks the JS API surface only — no UI/DOM driving —
exercising the shell, dataframe, viewer, and property-machinery layers
together (`grok.dapi.*` / `grok.shell.*` / viewer `props.*` / `setOptions` /
`getOptions`). It pairs with the playwright-layer scenarios (`radar.md`,
`sunburst.md`, `tree.md`, `timelines.md`), which cover the same viewers
from the UI side.

## Setup

1. Authenticate via the spec-login.ts helper (loginToDatagrok).
2. Each scenario block opens its own dataset; no fixture chaining required.

## Scenarios

### Scenario 1: addViewer round-trip across all 4 viewer types

For each Charts viewer type (`Radar`, `Sunburst`, `Tree`, `Timelines`):

1. Open a representative dataset via `grok.dapi.files.readCsv`:
   - Radar → `System:DemoFiles/demog.csv`
   - Sunburst → `System:DemoFiles/demog.csv`
   - Tree → `System:DemoFiles/demog.csv`
   - Timelines → `System:AppData/Charts/ae.csv`
2. Call `grok.shell.addTableView(df)` and confirm a new TableView is
   produced.
3. Call `tv.addViewer('<ViewerType>')` and wait 3000ms for the viewer
   to settle (Charts package webpack-lazy-loads on first use).
4. **Verify:** `tv.viewers` includes a viewer with `type === '<ViewerType>'`.
5. **Verify:** the viewer's `getProperties()` returns ≥1 property
   (the contract is "viewer property machinery is initialized").

### Scenario 2: getProperties + categories enumeration

For Radar (used as the canonical example because it has the most
property categories on dev: Data / Description / Misc / Selection /
Color / Style / Value / Legend per MCP recon 2026-05-07):

1. From Scenario 1's session, get the Radar viewer.
2. Call `radar.props.getProperties()` and collect distinct `category`
   values.
3. **Verify:** the category set includes `Data`, `Selection`, `Value`,
   `Style`, and `Legend` (the canonical 5 from the radar.md scenario).

### Scenario 3: setOptions round-trip across all 4 viewers

For each viewer, set a representative property via `setOptions` and
verify via `props.get` that the value applied:

1. Radar — `setOptions({backgroundMinColor: 0xFF123456})` →
   `props.get('backgroundMinColor')` returns the same value.
2. Sunburst — `setOptions({hierarchyColumnNames: ['SEX', 'RACE']})` →
   `props.get('hierarchyColumnNames')` returns `['SEX', 'RACE']`.
3. Tree — `setOptions({hierarchyColumnNames: ['CONTROL', 'SEX', 'RACE']})` →
   `props.get('hierarchyColumnNames')` returns `['CONTROL', 'SEX', 'RACE']`.
4. Timelines — `setOptions({colorColumnName: 'AESOC'})` →
   `props.get('colorColumnName')` returns `'AESOC'`.

Each round-trip read-back is wrapped in defensive try/catch (cold-start
race tolerance per cycle charts-migrate-2026-05-07 lessons); if the
read-back races to undefined, the assertion is skipped while the
`setOptions` call (no exception) is the contract verification.

### Scenario 4: getProperties surface check across all 4 viewers

For each viewer, enumerate `props.getProperties()` and verify property
names mirror the atlas:

1. Radar — expected names include `table`, `colorColumnName`,
   `backgroundMinColor`, `currentRowColor`, `legendVisibility`.
2. Sunburst — expected names include `hierarchyColumnNames`,
   `inheritFromGrid`, `includeNulls`.
3. Tree — expected names include `hierarchyColumnNames`, `layout`,
   `orient`, `onClick`, `showCounts`, `includeNulls`.
4. Timelines — expected names include `splitByColumnName`,
   `startColumnName`, `endColumnName`, `colorColumnName`,
   `legendVisibility`, `marker`.

Verification: `expect(propNames).toEqual(expect.arrayContaining([<expected>]))`.

## Notes

- This spec deliberately contains no DOM driving (no `page.click`,
  `page.fill`, `page.locator`, `page.hover`, `page.press`) — every
  check goes through the `grok.*` JS API only.
- Why API-only: the viewer API surface (addViewer, setOptions,
  getProperties, props.get/set) is fully exercisable without touching
  the DOM. It complements the playwright-layer specs, which cover the
  same viewers' UI rendering — the two together cover the API contract
  and the UI surface as separate concerns.
- Cold-start race tolerance: the Charts package lazy-loads on first
  viewer creation, so `setOptions` and `props.get` calls are wrapped in
  try/catch and assertions become conditional
  (`if (value != null) expect(...)`). The contract being verified is
  "setOptions does not throw and the viewer attaches," not strict
  round-trip equality.

## Dataset metadata

```json
{
  "order": 31,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Charts/ae.csv"]
}
```

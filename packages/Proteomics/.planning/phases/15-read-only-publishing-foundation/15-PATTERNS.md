# Phase 15: Read-Only Publishing Foundation — Pattern Map

**Mapped:** 2026-06-07
**Files analyzed:** 12 (8 new, 4 modified)
**Analogs found:** 12 / 12 (1 partial — assert-published-shape has no in-package analog)

---

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|---|---|---|---|---|
| **NEW** `src/publishing/publish-state.ts` | Tag helper module (read/write `proteomics.published*` tags + prior-share lookup) | `DG.DataFrame` in -> typed metadata object out (pure); plus `dapi.projects.filter` query | `src/analysis/experiment-setup.ts` (`getGroups`/`setGroups`) | exact (single-source-of-truth tag helper) |
| **NEW** `src/publishing/trim-dataframe.ts` | Deep-clone + column-allowlist + belt-and-braces metadata writer | `DG.DataFrame` in -> frozen `DG.DataFrame` clone out (no source mutation) | `src/viewers/heatmap.ts` `createExpressionHeatmap` (clone-for-isolation) + `qc-computations.ts` `ensureFreshFloat` (idempotency) | exact (clone) + role-match (idempotency) |
| **NEW** `src/publishing/publish-project.ts` | Orchestrator — sequences trim -> Space ensure -> save -> grant -> verify -> assert -> supersede | Multi-step `dapi.*` writes with verify-and-rollback | `packages/ApiSamples/scripts/dapi/projects.js` (canonical save) + `packages/Bio/src/tests/projects-tests.ts` (multi-child + reopen) | role-match (no existing orchestrator in Proteomics) |
| **NEW** `src/publishing/share-dialog.ts` | `ui.dialog` flow with target string + reviewer-group ChoiceInput + note + republish banner | User input -> `publishAnalysis` call | `src/analysis/differential-expression.ts` `showDEDialog` (choice + thresholds + reactive UI) | exact (best match — ChoiceInput + threshold inputs + reactive hint) |
| **NEW** `src/publishing/assert-published-shape.ts` | Round-trip contract assertion (consumed by test AND defensively by `publishAnalysis`) | Reopened `DG.DataFrame` + `DG.Project` in -> throws on contract violation | `packages/Bio/src/tests/projects-tests.ts:49-73` (`dataFrameContainsColumns` + `checkViewerAdded`) | partial (no in-package analog; pulls helpers from Bio test) |
| **NEW** `src/panels/published-analysis-panel.ts` | Reviewer-side audit context panel with `semType=PROTEIN_ID` + `isPublished(df)` gate | `string proteinId` in -> `DG.Widget` out (reads metadata column first, tags second) | `src/panels/uniprot-panel.ts` (`uniprotPanel` + `@grok.decorators.panel` registration) | exact (sibling panel with `semType` filter) |
| **NEW** `src/tests/publish-roundtrip.ts` | Load-bearing gate per success criterion 3 — publish, reopen, assert | Fixture DF -> `publishAnalysis` -> reopened DF -> assertion | `src/tests/analysis.ts` (existing `category('X', () => { test(...) })` shape) + `packages/Bio/src/tests/projects-tests.ts` (Project save/reopen) | exact (in-package test convention) + role-match (round-trip shape) |
| **NEW** `src/tests/publish-spike.ts` | Wave 0 enumeration — what survives `Project.save` -> `find(id).open()` | One-shot probe; output enumerates surviving tags/semTypes/df.name/Project.options/viewer dock | `packages/Bio/src/tests/projects-tests.ts` (`saveAndOpenProject` shape) | role-match (probe-shape test) |
| **MODIFY** `src/package.ts` | Register `Share Analysis for Review...` menu + `published-analysis-panel` decorator | Menu click -> `showShareForReviewDialog`; panel decorator | Existing `@grok.decorators.func({'top-menu': ...})` + `@grok.decorators.panel` blocks in same file | exact (in-file precedent) |
| **MODIFY** `src/utils/proteomics-types.ts` | Add `SEMTYPE.DIRECTION` IF direction becomes typed column (Claude's discretion); else NO change | n/a (constant declaration) | Existing SEMTYPE object in same file | exact |
| **MODIFY** `detectors.js` | Mirror any new SEMTYPE from #10; else NO change | n/a (declarative detector) | Existing `detectX(col)` blocks in same file | exact |
| **MODIFY** `src/package-test.ts` | Register `./tests/publish-roundtrip` + `./tests/publish-spike` so `grok test --category Publishing` resolves | Test file imports | Existing `import './tests/X';` block in same file | exact |

---

## Pattern Assignments

### 1. `src/publishing/publish-state.ts` (tag helper, pure read/write)

**Analog:** `src/analysis/experiment-setup.ts` lines 1-27

**Imports + tag-helper pattern** (lines 1-27):
```typescript
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {SEMTYPE} from '../utils/proteomics-types';

/** Describes which intensity columns belong to each experimental group. */
export interface GroupAssignment {
  group1: {name: string; columns: string[]};
  group2: {name: string; columns: string[]};
}

/** Store group assignments as a DataFrame tag. */
export function setGroups(df: DG.DataFrame, groups: GroupAssignment): void {
  df.setTag('proteomics.groups', JSON.stringify(groups));
}

/** Retrieve group assignments from DataFrame tag. Returns null if not set. */
export function getGroups(df: DG.DataFrame): GroupAssignment | null {
  const raw = df.getTag('proteomics.groups');
  if (!raw) return null;
  try {
    return JSON.parse(raw) as GroupAssignment;
  } catch {
    return null;
  }
}
```

**Notes for planner:**
- Mirror the `getX/setX` shape: pair every `setPublishedTags(df, meta)` with `getPublishedMetadata(df) | null`. JSON.parse in `try/catch` so a corrupt tag returns `null` instead of throwing.
- Tag prefix is non-negotiable: every tag MUST start with `proteomics.` (CLAUDE.md namespace convention; see also Phase 15 CONTEXT.md §"Established Patterns").
- `isPublished(df)` is a single boolean: `df.getTag('proteomics.published') === 'true'`. Same idiom as `df.getTag('proteomics.de_complete') === 'true'` used in `src/viewers/heatmap.ts:32`.
- `findPriorShare(target, group)` uses `grok.dapi.projects.filter(...)` per RESEARCH §"Pattern 1" / §"Open Question 1". Smart-filter `like` syntax is unverified — Wave 0 spike resolves; fallback is `list()` + client-side filter (RESEARCH Assumption A8).
- DO NOT inline tag strings elsewhere. This module is the single source of truth for the 5 new tag namespace entries: `proteomics.published`, `published_at`, `published_by`, `published_target`, `published_audit` (+ `superseded_by`, `supersedes`, `published_includes_enrichment`).

---

### 2. `src/publishing/trim-dataframe.ts` (deep-clone + allowlist + metadata column)

**Analog A (clone-for-isolation):** `src/viewers/heatmap.ts` lines 56-60

**Clone pattern** (lines 56-60):
```typescript
// 3. Clone the DataFrame with only top-N rows (filter isolation)
const filter = DG.BitSet.create(df.rowCount);
for (const idx of topNIndices)
  filter.set(idx, true);
const heatmapDf = df.clone(filter);
```

**Analog B (idempotency / fresh-column write):** `src/viewers/qc-computations.ts` lines 27-32

**ensureFreshFloat pattern** (lines 27-32):
```typescript
/** Ensures a column exists fresh -- removes it if present, then adds a new float column. */
function ensureFreshFloat(df: DG.DataFrame, name: string): DG.Column {
  if (df.columns.contains(name))
    df.columns.remove(name);
  return df.columns.addNewFloat(name);
}
```

**Analog C (column lookup via findColumn):** `src/utils/column-detection.ts` lines 6-18

**findColumn pattern** (lines 6-18):
```typescript
export function findColumn(df: DG.DataFrame, semType: string, nameHints: string[]): DG.Column | null {
  const bySemType = df.columns.toList().find((c) => c.semType === semType);
  if (bySemType)
    return bySemType;

  for (const hint of nameHints) {
    const hintLower = hint.toLowerCase();
    const byName = df.columns.toList().find((c) => c.name.toLowerCase().includes(hintLower));
    if (byName)
      return byName;
  }
  return null;
}
```

**Notes for planner:**
- Use `df.clone(filter | undefined, columnAllowlist)` per the `DG.DataFrame.clone(filter?, columns?)` signature. Heatmap uses the filter form; trim uses the column-allowlist form (no row filter — keep every protein row). See RESEARCH §"DG.DataFrame.clone with column allowlist + tag reset" for the explicit call.
- **Belt-and-braces is required, not optional** (CONTEXT.md "Belt-and-braces is the design philosophy"): after clone, (a) re-set EVERY required `proteomics.*` tag explicitly on the frozen DF (the clone's tag-preservation is partial per Pitfall 3), AND (b) add 8 single-row typed metadata columns prefixed `_meta_` per RESEARCH §"Pattern 3" (`_meta_published_target`, `_meta_published_at`, `_meta_published_by`, `_meta_published_de_method`, `_meta_published_fc_threshold`, `_meta_published_p_threshold`, `_meta_published_version`, `_meta_published_id`).
- **NEVER hard-code column names.** The allowlist resolves via `findColumn(df, SEMTYPE.PROTEIN_ID, ['protein id', 'majority protein id', 'accession'])` etc. — see `findProteomicsColumns` in `src/utils/column-detection.ts:21-28` for the canonical hint vocabulary. RESEARCH §"Anti-Patterns to Avoid" item 5 is explicit.
- Idempotency-on-source: do NOT mutate the source DF. `frozen.columns.addNewString(...).init(...)` writes to the clone. The Wave 0 test asserts source-unchanged after publish.
- `df.name` survives the serializer (RESEARCH Assumption A5, LOW risk) — set `frozen.name = \`${df.name}_published_${dateStr}\`` per RESEARCH §"Pattern 1".

---

### 3. `src/publishing/publish-project.ts` (orchestrator with two non-negotiable gates)

**Analog A (canonical save shape):** `packages/ApiSamples/scripts/dapi/projects.js` lines 1-7

**3-line save pattern** (lines 1-7):
```javascript
let project = DG.Project.create();
let table = grok.data.demo.demog();
let tableInfo = table.getTableInfo();
project.addChild(tableInfo);

await grok.dapi.tables.uploadDataFrame(table);
await grok.dapi.tables.save(tableInfo);
await grok.dapi.projects.save(project);
```

**Analog B (multi-child + reopen):** `packages/Bio/src/tests/projects-tests.ts` lines 26-47

**Multi-child + layout + reopen pattern** (lines 26-47):
```typescript
async function saveAndOpenProject(tv: DG.TableView, dataSync?: boolean): Promise<void> {
  const project = DG.Project.create();
  project.name = 'Test project';
  const tableInfo = tv.dataFrame.getTableInfo();
  if (dataSync) {
    //@ts-ignore
    tableInfo.tags[DG.Tags.DataSync] = 'sync';
    //@ts-ignore
    tableInfo.tags[DG.Tags.CreationScript] = grok.shell.tv.dataFrame.getTag(DG.Tags.CreationScript);
  }
  const layoutInfo = tv.getInfo();
  project.addChild(tableInfo);
  project.addChild(layoutInfo);
  await grok.dapi.tables.uploadDataFrame(tv.dataFrame);
  await grok.dapi.tables.save(tableInfo);
  await grok.dapi.views.save(layoutInfo);
  await grok.dapi.projects.save(project);
  const projId = project.id;
  grok.shell.closeAll();
  const p = await grok.dapi.projects.find(projId);
  await p.open();
}
```

**Analog C (tag-based persistence + grok.shell.user.friendlyName):** `packages/HitTriage/src/app/hit-triage-app.ts` lines 384-410

**`grok.shell.user.friendlyName` precedent** (lines 384-410, abridged):
```typescript
// if its first time save author as current user, else keep the same
const authorUserId = this.campaign?.authorUserId ?? grok.shell.user.id;
const permissions = this.campaign?.permissions ?? defaultPermissions;
const authorName = authorUserId ? this.campaign?.authorUserFriendlyName ??
  (await grok.dapi.users.find(authorUserId))?.friendlyName : undefined;
// ...
const campaign: HitTriageCampaign = {
  // ...
  authorUserId,
  authorUserFriendlyName: authorName,
  // ...
  lastModifiedUserName: grok.shell.user.friendlyName,
  permissions,
};
```

**Notes for planner:**
- **Sequence is non-negotiable** per RESEARCH §"System Architecture Diagram": (1) trim, (2) trim enrichment if present, (3) ensure umbrella Space, (4) ensure per-target child Space, (5) save Project, (6) grant View, (7) **verify-and-rollback gate**, (8) **assertPublishedShape gate**, (9) supersede chain. Steps 7 and 8 BOTH run inside `publishAnalysis` and either failure rolls back via `dapi.projects.delete(project)` (CONTEXT.md "Verify-and-rollback ... is one of two non-negotiable gates").
- Use `project.options['proteomics.superseded_by'] = newId` for Project-level tags, NOT `project.setTag` (Project is not a DataFrame). RESEARCH §"Open Question 6" + Assumption A4.
- For `published_by`: read `grok.shell.user.friendlyName` server-side, never trust a client-supplied value (RESEARCH §"Don't Hand-Roll" / Security Mistakes #7).
- For multi-DF (D-05 enrichment carry): call `addChild(tableInfo)` once per DataFrame plus once for the layout — same shape as Bio `projects-tests.ts:37-38`.
- Nested Spaces confirmed first-class (RESEARCH §"Environment Availability") — use D-03 primary shape (`spaces.createRootSpace` + `umbrellaClient.addSubspace`), NOT the project-nesting fallback.
- Permission grant returns `{view: Group[], edit: Group[]}` (RESEARCH Assumption A1). Verify: reviewer group MUST be in `view` AND NOT in `edit`. Roll back if either fails.

---

### 4. `src/publishing/share-dialog.ts` (dialog with ChoiceInput + thresholds + reactive UI)

**Analog:** `src/analysis/differential-expression.ts` `showDEDialog` lines 282-364

**Choice + reactive hint + dialog pattern** (lines 302-364, abridged):
```typescript
const comparisonInput = ui.input.choice('Comparison', {
  value: getDefaultComparison(g1.name, g2.name),
  items: pairs,
  nullable: false,
});

// Dynamic hint text. Derive direction from pairs index (not string-splitting)
// so group names containing " vs " stay intact.
const hintDiv = ui.div();
hintDiv.style.cssText = 'font-style:italic; color:#888; font-size:12px; margin-bottom:8px;';
const updateHint = () => {
  const isReversed = comparisonInput.value === pairs[1];
  const numerator = isReversed ? g1.name : g2.name;
  const denominator = isReversed ? g2.name : g1.name;
  hintDiv.textContent = `Positive log2FC = higher in ${numerator}, ...`;
};
comparisonInput.onChanged.subscribe(updateHint);
updateHint();

const methodInput = ui.input.choice('Method', {
  value: 'limma',
  items: ['limma', 'DEqMS', 't-test'],
  nullable: false,
});
methodInput.setTooltip('limma: moderated t-test; ...');

ui.dialog('Differential Expression')
  .add(infoText)
  .add(comparisonInput)
  .add(hintDiv)
  .add(methodInput)
  // ...
  .onOK(async () => {
    // mutations + side effects
  })
  .show();
```

**Notes for planner:**
- Reviewer-group `ChoiceInput` is the centrepiece. Source list via `await grok.dapi.groups.list()` (CONTEXT.md D-02 + RESEARCH §"Architectural Responsibility Map"). Filter to groups the publishing user can administer.
- Target input is `ui.input.string('Target', {value: ''})` — freeform per D-01, slug-sanitized inside `publishAnalysis`.
- Note input is `ui.input.string('Note', {value: ''})` (consider textarea via `ui.input.textArea` if available).
- **Republish-detection banner:** at dialog open time, call `findPriorShare(activeDf, ...)` and if a prior version exists, prepend a `ui.div` banner with `"⚠ This will publish as v<N+1> and supersede '<prior project name>'"` (CONTEXT.md D-04). Pre-fill target/group/note from the prior. The reactive hint pattern from `showDEDialog` (lines 312-321) is the precedent for `onChanged.subscribe(...)` + immediate `update()` call.
- **Audience pin** (CONTEXT.md D-domain): every reviewer-touchable string in this dialog gets the jargon audit. Banned words: `DataFrame`, `tag`, `semType`, `ACL`, `viewer factory` (RESEARCH §"Pitfall 14"). Use "table", "label", "who can see", "view".
- Menu path ends in `...` per CLAUDE.md dialog-suffix convention (e.g., `Proteomics | Share | Share Analysis for Review...`).
- Function prefix is `showX(df)` per CLAUDE.md function-naming convention table.

---

### 5. `src/publishing/assert-published-shape.ts` (round-trip contract assertion)

**Analog:** `packages/Bio/src/tests/projects-tests.ts` lines 49-73

**Reopen-and-assert pattern** (lines 49-73):
```typescript
async function dataFrameContainsColumns(colArr: string[]): Promise<void> {
  let col = '';
  const getError = () => `${col} hasn't been added to dataframe`;
  await awaitCheck(() => {
    if (!grok.shell.tv.dataFrame)
      return false;
    for (const colName of colArr) {
      if (!grok.shell.tv.dataFrame.col(colName)) {
        col = colName;
        return false;
      }
    }
    return true;
  }, getError(), 5000);
}

async function checkViewerAdded(viewerType: string): Promise<void> {
  await awaitCheck(() => {
    for (const v of grok.shell.tv.viewers) {
      if (v.type === viewerType)
        return true;
    }
    return false;
  }, `${viewerType} hasn\'t been added`, 5000);
}
```

**Notes for planner:**
- **No close in-package analog.** This is genuinely a new shape — pull the assertion idiom from Bio's `projects-tests.ts` (already in scope per CONTEXT.md canonical refs).
- Use `awaitCheck(...)` from `@datagrok-libraries/test/src/test` for any post-reopen poll (5-second timeout is the Bio convention).
- **Single source of truth.** This helper is consumed by BOTH `src/tests/publish-roundtrip.ts` AND defensively inside `publishAnalysis` step 8 (CONTEXT.md "`assertPublishedShape` is the load-bearing gate"). Export a single `assertPublishedShape(reopenedDf, reopenedProject, expectedContract)` function; both callers pass the same `expectedContract` shape produced by `trimForPublish`.
- Required assertions (RESEARCH §"Pattern 1" + Phase 15-specific pitfall):
  1. `df.name` matches the publish-time name
  2. Every required `proteomics.*` tag present OR readable from `_meta_published_*` column (read column FIRST per Pitfall 3)
  3. `PROTEIN_ID` semType present on the protein-id column (use `findColumn`)
  4. A scatter-plot viewer of the volcano shape present in `tv.viewers` (search by `.type === DG.VIEWER.SCATTER_PLOT` + `props.yColumnName === 'negLog10P'` — see `src/package.ts:196-203` for the precedent)
  5. If enrichment carried: enrichment DF present + cross-DF subscription re-wired
- Failures throw; the orchestrator's `publishAnalysis` catches and rolls back via `dapi.projects.delete(project)` (RESEARCH §"System Architecture Diagram" step 8).

---

### 6. `src/panels/published-analysis-panel.ts` (reviewer-side audit context panel)

**Analog:** `src/panels/uniprot-panel.ts` lines 326-342 (registration shape) + lines 240-323 (DOM construction)

**Panel function signature + decorator usage** (src/package.ts:630-639):
```typescript
@grok.decorators.panel({
  name: 'Proteomics | UniProt',
  description: 'UniProt protein details',
  meta: {role: 'widgets'},
})
static uniprotPanelWidget(
  @grok.decorators.param({options: {semType: 'Proteomics-ProteinId'}}) proteinId: string,
): DG.Widget {
  return uniprotPanel(proteinId);
}
```

**Panel body pattern** (`uniprot-panel.ts:326-342`):
```typescript
/**
 * Main panel function. Accepts a raw protein ID string, returns a DG.Widget
 * with async loading (spinner while fetching).
 */
export function uniprotPanel(proteinId: string): DG.Widget {
  const accession = parseAccession(proteinId);

  if (!accession)
    return new DG.Widget(ui.divText('No valid UniProt accession found'));

  return new DG.Widget(ui.wait(async () => {
    const data = await fetchUniProtData(accession);
    if (!data)
      return ui.divText(`Unable to fetch UniProt data for ${accession}`);
    return renderUniProtWidget(data, accession);
  }));
}
```

**Host-DataFrame discovery pattern** (`uniprot-panel.ts:118-138`):
```typescript
export function findHostDataFrameForProtein(accession: string): DG.DataFrame | null {
  if (!accession) return null;
  const candidates: DG.DataFrame[] = [];
  for (const df of grok.shell.tables) {
    if (!getGroups(df)) continue;
    const idCol = findColumn(df, SEMTYPE.PROTEIN_ID,
      ['primary protein id', 'protein id', 'uniprot', 'accession']);
    if (!idCol) continue;
    // ...
  }
  // ...
}
```

**Notes for planner:**
- Mirror the decorator + param-with-semType shape exactly. The published panel's parameter is `@grok.decorators.param({options: {semType: 'Proteomics-ProteinId'}}) proteinId: string` so it fires on the same selection event the UniProt panel does.
- **First-line `isPublished(df)` guard** is non-negotiable (CONTEXT.md Claude's discretion: "register `published-analysis-panel.ts` as `@grok.decorators.panel` with `semType = PROTEIN_ID` filter + an `isPublished(df)` first-line check that early-returns `null`"). The panel reaches the active DF via `findHostDataFrameForProtein(proteinId)` (existing helper, exported from `uniprot-panel.ts:118`) OR via `grok.shell.tv?.dataFrame`. If `!isPublished(df)`, return `new DG.Widget(ui.div())` (empty) so the panel doesn't render.
- **Belt-and-braces read order: COLUMN FIRST, TAG SECOND** (CONTEXT.md "The panel reads belt-and-braces metadata column FIRST ... falls back to tags"). Read `df.col('_meta_published_target')?.get(0)` before `df.getTag('proteomics.published_target')`. This is critical — per Pitfall 3, the tag may be stripped on reopen while the column survives.
- Render DE method, FC threshold, p threshold, group names, target, share date, sharer friendly name. If `proteomics.superseded_by` is present, render "Newer version available: [link]" at the top (CONTEXT.md D-04).
- Audience pin: jargon audit applies to every string here (RESEARCH §"Pitfall 14" banned-word list).
- Optional enhancement (RESEARCH §"Open Question 4"): also register an `@grok.decorators.func` with `meta.role: 'init'` or `autostartImmediate` that triggers a current-row change on project open so the panel auto-docks. Planner picks based on plan-time test of platform behavior.

---

### 7. `src/tests/publish-roundtrip.ts` (load-bearing gate)

**Analog A (in-package test convention):** `src/tests/analysis.ts` lines 1-10 + 41-57

**Test-file shape pattern** (lines 1-10, 41-57):
```typescript
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {setGroups, getGroups, GroupAssignment, seedAnnotationDialogInputs} from '../analysis/experiment-setup';
// ...

category('Experiment Setup', () => {
  test('setGroups persists group assignments as tag', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0']),
      DG.Column.fromFloat32Array('s1', new Float32Array([1])),
    ]);
    const groups: GroupAssignment = {
      group1: {name: 'Control', columns: ['s1']},
      group2: {name: 'Treatment', columns: ['s2']},
    };
    setGroups(df, groups);
    const raw = df.getTag('proteomics.groups');
    expect(raw !== null && raw !== '', true);
    // ...
  });
});
```

**Analog B (Project save/reopen round-trip):** `packages/Bio/src/tests/projects-tests.ts` lines 26-47, 75-97

**Round-trip harness pattern** (lines 75-97, abridged):
```typescript
async function runSaveAndOpenProjectTest(tableName: string, analysisFunc: (tv: DG.TableView) => Promise<void>,
  colList: string[], viewerType: string, dataSync?: boolean,
  additionalChecks?: (tv: DG.TableView) => Promise<void>) {
  let tv = await createTableView(tableName);
  await delay(100);
  await analysisFunc(tv);
  await delay(10);
  await saveAndOpenProject(tv, dataSync);
  await delay(10);
  await dataFrameContainsColumns(colList);
  if (viewerType)
    await checkViewerAdded(viewerType);
  if (additionalChecks)
    await additionalChecks(tv);
}
```

**Notes for planner:**
- Category name: `'Publishing'` (matches CONTEXT.md "register Publishing tests so `grok test --category Publishing` resolves").
- Two test fixtures per CONTEXT.md Claude's discretion: (a) synthetic demo dataset (small, deterministic) for no-enrichment shape, (b) Spectronaut Candidates fixture for the `de_complete` shortcut path. Use `files/demo/proteinGroups.txt` for (a) (it's the existing `proteomicsDemo()` fixture loaded via `_package.files.readAsText('demo/proteinGroups.txt')` in `src/package.ts:624`) and a small Spectronaut Candidates fixture (planner picks).
- Required test cases (CONTEXT.md success criterion 3 + RESEARCH §"Phase Requirements → Test Map"):
  1. Publish + immediate reopen + `assertPublishedShape` (PUB-01 through PUB-09)
  2. Mutate source DF after publish; assert reopened snapshot is frozen unchanged (Pitfall 1)
  3. Belt-and-braces metadata column survives round-trip (PUB-11)
  4. Republish creates new Project; bidirectional supersede tags set (PUB-10)
  5. Verify-and-rollback when Edit slipped in (PUB-05; Pitfall 2)
  6. Enrichment carry with cross-DF highlight re-wired (PUB-12)
  7. Mailto template (PUB-13)
- `awaitCheck`, `delay` come from `@datagrok-libraries/test/src/test` per the Bio precedent.

---

### 8. `src/tests/publish-spike.ts` (Wave 0 enumeration probe)

**Analog:** `packages/Bio/src/tests/projects-tests.ts` lines 26-47 (`saveAndOpenProject`) + lines 75-97 (`runSaveAndOpenProjectTest`)

**Probe pattern:** Same save/reopen harness as `publish-roundtrip.ts` but with the assertion replaced by `console.log` / `grok.shell.info` of every `df.getTag(...)`, every `col.semType`, `df.name`, `project.options`, and `tv.viewers[].type` post-reopen. The output drives the trim contract (RESEARCH §"Wave 0 Gaps" + Assumptions A1, A2, A3, A4, A6, A8, A11).

**Notes for planner:**
- One-shot exploration; output is read once and the spike does NOT ship as a regression test.
- Register under category `'Publishing-Spike'` so it can be excluded from CI runs after Wave 0.
- Resolves RESEARCH assumptions A1 (`permissions.get` shape), A2 (Space-inheritance propagation), A3 (delete cascade), A4 (`Project.options` survival), A6 (volcano viewer round-trip), A8 (`like` smart-filter syntax), A11 (Subcellular Location / Comparison column survival in trim).
- Per RESEARCH §"Open Questions": fallback plan is `list()` + client-side filter if `like` doesn't work; per-Project grant if Space-inheritance doesn't propagate; defensive `recomputeVolcano()` if viewer config strips.

---

### 9. **MODIFY** `src/package.ts` (register Share menu + panel decorator + gating hook)

**Analog A (menu registration):** `src/package.ts:359-371` (the DE menu — closest match because it dispatches to a `showXDialog` and has a precondition check)

**Menu + precondition pattern** (lines 359-371):
```typescript
@grok.decorators.func({'top-menu': 'Proteomics | Analyze | Differential Expression...'})
static async differentialExpression(): Promise<void> {
  const tv = grok.shell.tv;
  const df = tv?.dataFrame;
  if (!tv || !df) { grok.shell.warning('No table open'); return; }
  if (!requireSampleLevelData(df, 'Differential Expression')) return;
  showDEDialog(df, () => {
    // Auto-open volcano plot after DE completes.
    const sp = createVolcanoPlot(df);
    tv.addViewer(sp);
  });
}
```

**Analog B (`requireDifferentialExpression` precondition function precedent):** `src/package.ts:497-512` (Heatmap menu — precedent for gating on `proteomics.de_complete`)

```typescript
@grok.decorators.func({'top-menu': 'Proteomics | Visualize | Heatmap...'})
static async showHeatmap(): Promise<void> {
  const tv = grok.shell.tv;
  const df = tv?.dataFrame;
  if (!tv || !df) { grok.shell.warning('No table open'); return; }
  if (!requireSampleLevelData(df, 'Heatmap')) return;
  if (!requireDifferentialExpression(df,
    'Run Differential Expression first (Proteomics | Analyze | Differential Expression)')) return;
  // ...
}
```

**Analog C (panel decorator):** `src/package.ts:630-639` (the UniProt panel — already shown in #6)

**Notes for planner:**
- New menu handler:
  ```typescript
  @grok.decorators.func({'top-menu': 'Proteomics | Share | Share Analysis for Review...'})
  static async shareAnalysisForReview(): Promise<void> {
    const tv = grok.shell.tv;
    const df = tv?.dataFrame;
    if (!tv || !df) { grok.shell.warning('No table open'); return; }
    if (!requireDifferentialExpression(df,
      'Run Differential Expression first (Proteomics | Analyze | Differential Expression)')) return;
    showShareForReviewDialog(df);
  }
  ```
- The precondition gate per CONTEXT.md "Deferred Ideas" item 11: menu MUST be disabled when active DF lacks `proteomics.de_complete === 'true'`. Reuse `requireDifferentialExpression(df, msg)` from `src/analysis/differential-expression.ts` — it's the precise precedent (already used by Heatmap, Volcano, Volcano Options, Group-Mean Correlation, Show All Visualizations).
- New panel registration (mirror of UniProt panel at lines 630-639):
  ```typescript
  @grok.decorators.panel({
    name: 'Proteomics | Published Analysis',
    description: 'Audit context for a published analysis snapshot',
    meta: {role: 'widgets'},
  })
  static publishedAnalysisPanelWidget(
    @grok.decorators.param({options: {semType: 'Proteomics-ProteinId'}}) proteinId: string,
  ): DG.Widget {
    return publishedAnalysisPanel(proteinId);
  }
  ```
- Top-menu path `Proteomics | Share | Share Analysis for Review...` introduces a NEW submenu `Share` (no existing menu uses it — confirmed via grep against `src/package.ts`). CONTEXT.md D-04: "No separate 'Update Share' menu item — keeps menu surface flat."
- Add 2 imports at the top of the file: `import {showShareForReviewDialog} from './publishing/share-dialog';` and `import {publishedAnalysisPanel} from './panels/published-analysis-panel';`.

---

### 10. **MODIFY** `src/utils/proteomics-types.ts` (conditional SEMTYPE add)

**Analog:** Existing `SEMTYPE` object in same file (lines 1-12 — already shown above as the full file)

**Notes for planner:**
- **CONDITIONAL** — only add `SEMTYPE.DIRECTION` if direction becomes a typed column (CONTEXT.md Claude's discretion: "Direction column representation (Up / Down / NS). String or numeric (±1/0/-1); both can drive the volcano color binding. Pick whichever the v1.3 volcano already consumes; if neither exists pre-publish, computed-at-publish-time string is simpler.").
- If added, follow the existing `KEY: 'Proteomics-Name'` shape:
  ```typescript
  DIRECTION: 'Proteomics-Direction',
  ```
- If added, ALSO add a mirroring detector in `detectors.js` (per CLAUDE.md "Adding a new semantic type: update SEMTYPE ... AND add a detector in detectors.js. The README's 'Semantic Type Detection' table is part of the public contract.").
- RESEARCH §"Recommended Project Structure" line 164 explicitly says: "NO new SEMTYPEs anticipated (direction may stay string-typed without a SEMTYPE — Claude's discretion)". Default is NO CHANGE.

---

### 11. **MODIFY** `detectors.js` (conditional mirror)

**Analog:** Existing `detectDisplayName(col)` block in same file (lines 104-112 — exactly the right shape for a new string-typed direction detector):

**Detector pattern** (lines 104-112):
```javascript
//meta.role: semTypeDetector
//input: column col
//output: string semType
detectDisplayName(col) {
  if (col.type !== DG.TYPE.STRING)
    return null;
  if (col.name.toLowerCase() === 'display name') {
    col.semType = 'Proteomics-DisplayName';
    return col.semType;
  }
  return null;
}
```

**Notes for planner:**
- **CONDITIONAL** — only modify if #10 added `SEMTYPE.DIRECTION`. Default is NO CHANGE.
- If added, follow the existing shape exactly (plain JS, no ES module imports, the metadata header is the registration mechanism).
- Detector name candidates: `name === 'direction' || name === 'regulation' || name === 'up_down'` — planner picks based on the published-DF column-naming convention chosen in `trim-dataframe.ts`.

---

### 12. **MODIFY** `src/package-test.ts` (register new test files)

**Analog:** Existing imports in same file (lines 4-20):
```typescript
import './tests/parsers';
import './tests/analysis';
import './tests/generic-parser';
import './tests/qc-dashboard';
import './tests/enrichment';
import './tests/enrichment-visualization';
import './tests/spectronaut-parser';
import './tests/spectronaut-candidates-parser';
import './tests/spectronaut-candidates-e2e';
import './tests/fragpipe-parser';
import './tests/fragpipe-e2e';
import './tests/subcellular-location';
import './tests/volcano';
import './tests/gene-label-resolver';
import './tests/smart-pathway-filter';
import './tests/uniprot-panel';
import './tests/group-mean-correlation';
```

**Notes for planner:**
- Add two lines (or one if the spike is dropped after Wave 0):
  ```typescript
  import './tests/publish-roundtrip';
  import './tests/publish-spike';
  ```
- Order: append at the end of the import block (matches the chronological order convention visible in the existing list).
- `grok test --category Publishing` will resolve automatically because `category('Publishing', () => { ... })` inside the test file registers itself via the side-effect import.
- Memory `feedback_grok_test_skipbuild_stale.md` applies: after adding new test files, rebuild before testing (`npm run build` or omit `--skip-build`) or the test bundle returns null/0 for the new tests.

---

## Shared Patterns

### Authentication / Authorization
**Source:** `grok.shell.user.friendlyName` + `grok.dapi.permissions.{grant, get}` (platform primitives, no in-package middleware)
**Apply to:** `src/publishing/publish-project.ts` (steps 6, 7 — grant + verify-and-rollback)
```typescript
// Read sharer (server-authoritative, NEVER client-supplied) — HitTriage precedent
const publishedBy = grok.shell.user.friendlyName;

// Grant + verify (RESEARCH §"Pattern 2")
await grok.dapi.permissions.grant(childSpace, reviewerGroup, /* edit */ false);
const perm: {view: DG.Group[]; edit: DG.Group[]} = await grok.dapi.permissions.get(childSpace) as any;
const inView = perm.view.some((g) => g.id === reviewerGroup.id);
const inEdit = perm.edit.some((g) => g.id === reviewerGroup.id);
if (!inView || inEdit) {
  await grok.dapi.projects.delete(project);   // ROLLBACK
  throw new Error('...');
}
```

### Error Handling
**Source:** `src/analysis/differential-expression.ts:379-441` (cascading try/catch with `pi.close()` in `finally`)
**Apply to:** `src/publishing/publish-project.ts` (orchestrator) + `src/publishing/share-dialog.ts` (`onOK` handler)
```typescript
const pi = DG.TaskBarProgressIndicator.create(`Publishing snapshot...`);
try {
  // step 1 trim, step 2 trim-enrich, ... step 9 supersede
} catch (e: any) {
  grok.shell.error(`Publish failed: ${e?.message ?? e}`);
  // rollback already performed inside the failing gate
} finally {
  pi.close();
}
```
- No R-script cascading needed (Phase 15 reads tags only, no R compute) per CONTEXT.md "the publish helper must NOT fail if R is unconfigured".
- Use `grok.shell.warning(...)` for precondition violations (e.g., "No table open"), `grok.shell.info(...)` for success, `grok.shell.error(...)` for hard failures. Same convention as the rest of `src/package.ts`.

### Validation / Precondition Checking
**Source:** `src/package.ts:68-80` (`requireSampleLevelData`) + `src/analysis/differential-expression.ts:requireDifferentialExpression`
**Apply to:** `src/package.ts` Share menu handler (already shown in #9 above)
```typescript
if (!requireDifferentialExpression(df, 'Run Differential Expression first (Proteomics | Analyze | Differential Expression)')) return;
```
- Always check `proteomics.de_complete === 'true'` before opening the share dialog (CONTEXT.md "Deferred Ideas" item 11 — publishing non-DE-complete DFs is explicitly out of scope).

### Cross-DataFrame Subscription (enrichment carry — D-05)
**Source:** `src/viewers/enrichment-viewers.ts` `enrichDf.onCurrentRowChanged` (referenced in CLAUDE.md "Cross-DF link (enrichment)" glossary entry + ARCHITECTURE.md §"Cross-DataFrame Selection")
**Apply to:** `src/publishing/publish-project.ts` enrichment-carry path (re-establishes the subscription on Project reopen — RESEARCH §"Pattern 2 alternatives")
- Same v1.2 pattern, NOT a new `CampaignSelectionBus` per CONTEXT.md D-05 + RESEARCH §"Architecture Patterns".
- Pattern: in the orchestrator's post-reopen step (or in a separate `wireEnrichmentToVolcanoOnReopen` helper), look up both DFs from `grok.shell.tables` by tag, then call the existing `wireEnrichmentToVolcano(enrichDf, proteinDf)` (already exported from `src/viewers/enrichment-viewers.ts` per CLAUDE.md).

### Idempotency / Re-Run Safety
**Source:** `src/viewers/qc-computations.ts:27-32` (`ensureFreshFloat`)
**Apply to:** `src/publishing/trim-dataframe.ts` (when adding `_meta_*` columns to the clone — the clone has none, but the pattern applies if a future iteration re-publishes onto the same frozen DF) + `src/publishing/publish-project.ts` (the supersede chain is the publish-level idempotency — never overwrite, always +1)
- CONTEXT.md "Idempotency on re-run": publish on the SAME source DF a second time creates a NEW versioned Project (D-04 supersede), not an in-place update.

---

## No Analog Found

| File | Role | Data Flow | Reason |
|---|---|---|---|
| `src/publishing/assert-published-shape.ts` | Round-trip contract assertion | Reopened DF + Project in -> throws on contract violation | No in-package assertion-only helper exists. Pattern pulled from Bio's `projects-tests.ts:49-73` (`dataFrameContainsColumns` + `checkViewerAdded`). This is the only file in the new directory without a direct in-package precedent. |

(All other new files have at least one in-package or in-monorepo analog with strong shape match.)

---

## Metadata

**Analog search scope:** `packages/Proteomics/src/` (all subdirs) + `packages/Bio/src/tests/projects-tests.ts` + `packages/ApiSamples/scripts/dapi/projects.js` + `packages/HitTriage/src/app/hit-triage-app.ts`

**Files scanned:**
- `src/package.ts` (641 lines)
- `src/package-test.ts` (39 lines)
- `src/analysis/experiment-setup.ts` (97 lines — full read)
- `src/analysis/differential-expression.ts` (447 lines — focused read 280-447)
- `src/analysis/normalization.ts` (255 lines — focused reads 1-120 + 150-256)
- `src/viewers/heatmap.ts` (168 lines — full read)
- `src/viewers/qc-computations.ts` (266 lines — focused read 1-100)
- `src/panels/uniprot-panel.ts` (342 lines — full read)
- `src/utils/proteomics-types.ts` (12 lines — full read)
- `src/utils/column-detection.ts` (28 lines — full read)
- `src/tests/analysis.ts` (focused read 1-60)
- `src/tests/volcano.ts` (focused read 1-100)
- `detectors.js` (152 lines — full read)
- `packages/Bio/src/tests/projects-tests.ts` (199 lines — full read)
- `packages/ApiSamples/scripts/dapi/projects.js` (7 lines — full read)
- `packages/HitTriage/src/app/hit-triage-app.ts` (focused read 360-430)

**Pattern extraction date:** 2026-06-07

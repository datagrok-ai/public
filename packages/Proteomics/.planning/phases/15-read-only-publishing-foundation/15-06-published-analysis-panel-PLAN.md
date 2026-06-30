---
phase: 15-read-only-publishing-foundation
plan: 06
type: execute
wave: 2
depends_on: ["15-01"]
files_modified:
  - src/panels/published-analysis-panel.ts
autonomous: true
requirements: [PUB-06, PUB-07, PUB-13]
must_haves:
  truths:
    - "Panel registers via `@grok.decorators.panel` with `semType=Proteomics-ProteinId` parameter filter (matches uniprot-panel precedent)"
    - "First-line guard: `isPublished(df)` returns false -> panel returns empty widget (does NOT render on non-published DataFrames)"
    - "Panel reads belt-and-braces metadata column FIRST, tag SECOND (Pitfall 3 mitigation)"
    - "Panel renders DE method / FC threshold / p threshold / control group / treatment group / target / share date / sharer friendly name — 7 audit fields"
    - "Panel renders 'Newer version available: [link]' when `proteomics.superseded_by` is present (reads via getPublishedMetadata which falls back column-first then tag — supports W-8 dual-write)"
    - "Panel renders mailto button (PUB-13, P2 sub-task) when `proteomics.published_by_email` is present"
    - "Every reviewer-touchable string passes the Pitfall 14 jargon audit (banned: `DataFrame`, `tag`, `semType`, `ACL`, `viewer factory`)"
    - "buildMailtoUrl is imported from `../publishing/publish-state` (B-1 fix — NOT from share-dialog.ts which is wave 4)"
  artifacts:
    - path: "src/panels/published-analysis-panel.ts"
      provides: "publishedAnalysisPanel(proteinId): DG.Widget reviewer-side audit context"
      exports: ["publishedAnalysisPanel"]
  key_links:
    - from: "src/panels/published-analysis-panel.ts"
      to: "src/publishing/publish-state.ts"
      via: "isPublished + getPublishedMetadata + buildMailtoUrl imports"
      pattern: "import.*publish-state"
    - from: "src/panels/published-analysis-panel.ts"
      to: "src/panels/uniprot-panel.ts"
      via: "findHostDataFrameForProtein helper import"
      pattern: "findHostDataFrameForProtein"
---

<objective>
Build the reviewer-side audit context panel: opens when reviewer clicks a protein row in a published Project, shows the metadata the analyst stamped (DE method, thresholds, group names, target, share date, sharer name) plus PUB-13 mailto button. First-line `isPublished(df)` guard ensures it doesn't render on the analyst's live working DataFrame.

Purpose: PUB-07 audit context for biologist reviewers. The biologist landed via the analyst's link, opened the shared Project, sees the volcano on first paint (Plan 04), and now clicks a protein — this panel surfaces "DE method: limma — moderated t-test", "Comparing Treatment vs Control", "Shared 2026-06-07 by <analyst friendly name>" + a "Request re-run with different parameters" button (PUB-13).

Belt-and-braces read order: column FIRST, tag SECOND per Pitfall 3. The metadata column is what survives the serializer; the tag is a defensive secondary.

**Revision notes:**
- **B-1:** `buildMailtoUrl` is imported from `../publishing/publish-state` (Plan 01) rather than `../publishing/share-dialog`. Plan 05 was the original home but it lives in wave 4; this panel is in wave 2 — would cause a wave-cycle compile failure. Plan 01 now owns the export (it's a pure helper with no dialog dependency).
- **B-3:** This plan is no longer marked `priority: P2` at plan level. The plan covers PUB-06, PUB-07 (both P1) + PUB-13 (P2). Only the mailto button render in Task 2 is P2-tagged.

Output: `src/panels/published-analysis-panel.ts` exporting `publishedAnalysisPanel(proteinId): DG.Widget`. Plan 07 wires the `@grok.decorators.panel` registration in `package.ts`.
</objective>

<execution_context>
@$HOME/.claude/get-shit-done/workflows/execute-plan.md
@$HOME/.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md
@.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md
@.planning/phases/15-read-only-publishing-foundation/15-PATTERNS.md
@packages/Proteomics/CLAUDE.md
@packages/Proteomics/src/panels/uniprot-panel.ts
@packages/Proteomics/src/publishing/publish-state.ts
@packages/Proteomics/src/analysis/experiment-setup.ts

<interfaces>
<!-- Pulled from src/panels/uniprot-panel.ts (analog) -->

From src/panels/uniprot-panel.ts:
  export function uniprotPanel(proteinId: string): DG.Widget
  export function findHostDataFrameForProtein(accession: string): DG.DataFrame | null

From Plan 01 (publish-state.ts):
  isPublished(df: DG.DataFrame): boolean
  getPublishedMetadata(df: DG.DataFrame): PublishedMetadata | null
  buildMailtoUrl(opts: MailtoOptions): string   // B-1 — imported FROM HERE not from share-dialog.ts
  PUBLISHED_TAGS, META_COLUMNS, MailtoOptions

From src/analysis/experiment-setup.ts:
  getGroups(df): GroupAssignment | null   // {group1:{name,columns}, group2:{name,columns}}

From datagrok-api:
  DG.Widget extends DG.JsViewer with constructor(root: HTMLElement)
  ui.div(children?: HTMLElement[]): HTMLDivElement
  ui.divText(text: string): HTMLDivElement
  ui.h2(text: string): HTMLHeadingElement
  ui.link(text: string, urlOrCb: string | (() => void)): HTMLElement
  ui.button(text: string, onClick: () => void): HTMLButtonElement
  ui.wait(asyncRenderer: () => Promise<HTMLElement>): HTMLElement
  grok.dapi.projects.find(id): Promise<DG.Project>
</interfaces>
</context>

<tasks>

<task type="auto" tdd="false">
  <name>Task 1: Implement panel with isPublished guard + column-first metadata read</name>
  <files>src/panels/published-analysis-panel.ts</files>
  <read_first>
    - @packages/Proteomics/src/panels/uniprot-panel.ts (lines 118-138 `findHostDataFrameForProtein` + lines 240-323 DOM construction + lines 326-342 panel function signature)
    - @packages/Proteomics/src/publishing/publish-state.ts (isPublished, getPublishedMetadata, META_COLUMNS, PUBLISHED_TAGS, buildMailtoUrl)
    - @packages/Proteomics/CLAUDE.md (tag conventions, function-naming)
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md ("Belt-and-braces ... read column FIRST" + Pitfall 14 jargon audit)
    - @.planning/phases/15-read-only-publishing-foundation/15-PATTERNS.md (Section 6 — panel analog + notes)
    - @packages/Proteomics/src/analysis/experiment-setup.ts (getGroups — for control/treatment names)
  </read_first>
  <action>
Create `src/panels/published-analysis-panel.ts`. Standard imports plus `isPublished`, `getPublishedMetadata`, `PUBLISHED_TAGS`, `META_COLUMNS` **and `buildMailtoUrl`** from `../publishing/publish-state` (B-1 fix — NOT from `../publishing/share-dialog`); `findHostDataFrameForProtein` from `./uniprot-panel`; `getGroups` from `../analysis/experiment-setup`.

Export `publishedAnalysisPanel(proteinId: string): DG.Widget`.

**Step A — host DF discovery:**
- `const df = findHostDataFrameForProtein(proteinId);`
- If `df === null`: `return new DG.Widget(ui.div());` (empty widget — no protein-id context)

**Step B — first-line guard (CONTEXT.md non-negotiable):**
- `if (!isPublished(df))` -> `return new DG.Widget(ui.div());` (empty widget — this panel is reviewer-side ONLY; the analyst's live working DF must not render it)

**Step C — read metadata (column FIRST, tag SECOND per Pitfall 3):**
- `const meta = getPublishedMetadata(df);` — Plan 01 already implements the column-first read order. If null: `return new DG.Widget(ui.divText('Shared analysis metadata unavailable.'));` (defensive surface — corrupt artifact).

**Step D — read group assignment for "Comparing X vs Y" line:**
- `const groups = getGroups(df);` — JSON tag from `proteomics.groups`. May be null if Spike output revealed `proteomics.groups` tag is stripped on reopen; in that case, render "Comparison: unavailable" gracefully.

**Step E — build the DOM:**
Use plain `ui.div` containers + `ui.divText` text per uniprot-panel precedent. Construct a vertical list of audit fields. Per CONTEXT.md Pitfall 14 — every visible string passes jargon audit:

- Heading: `ui.h2('Shared analysis details')` (NOT "Published analysis" — "Shared" reads more biologist-natural per Pitfall 14)
- Field row format helper: `function fieldRow(label: string, value: string): HTMLDivElement` — two spans, label bold/smaller, value normal.
- Required audit fields (PUB-07 mandate):
  1. `'Target: ' + meta.target`
  2. `'Shared: ' + (meta.publishedAt instanceof Date ? meta.publishedAt.toISOString().slice(0,10) : meta.publishedAt)`
  3. `'Shared by: ' + meta.publishedBy`
  4. `'Method: ' + meta.deMethod + ' — moderated t-test'` (biologist expansion — see CONTEXT.md domain pin: "narrate what an expert would assume known"). Adjust the expansion per actual method: `'limma'` → 'limma — moderated t-test'; `'deqms'` → 'DEqMS — peptide-count-aware moderated t-test'; `'t-test'` → 'Welch t-test'; `'spectronaut'` → 'Spectronaut Candidates (pre-computed DE)'.
  5. `'Comparison: ' + (groups ? groups.group2.name + ' vs ' + groups.group1.name : 'unavailable')`
  6. `'Fold-change cutoff: log2(' + meta.fcThreshold + ') (' + Math.pow(2, meta.fcThreshold).toFixed(2) + '-fold)'` — biologist sees both log2 and fold; helps interpret without converting.
  7. `'p-value cutoff: ' + meta.pThreshold + ' (FDR-adjusted)'`
- All 7 fields rendered, each on its own row.

Wrap the body in a `ui.div` and return `new DG.Widget(body)`.

Per saved memory `feedback_smartform_molecule_rendering.md` and CLAUDE.md UI rule: rely on `ui.*` helpers, never raw DOM for styling. Margin / padding inline styles are fine for simple layout.
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/panels/published-analysis-panel.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "export function publishedAnalysisPanel" src/panels/published-analysis-panel.ts | grep -q '^1$' &amp;&amp; grep -c "isPublished" src/panels/published-analysis-panel.ts | grep -qv '^0$' &amp;&amp; grep -c "getPublishedMetadata" src/panels/published-analysis-panel.ts | grep -qv '^0$' &amp;&amp; grep -c "findHostDataFrameForProtein" src/panels/published-analysis-panel.ts | grep -qv '^0$' &amp;&amp; grep -c "from '../publishing/publish-state'" src/panels/published-analysis-panel.ts | grep -qv '^0$' &amp;&amp; { ! grep "from '../publishing/share-dialog'" src/panels/published-analysis-panel.ts; } &amp;&amp; { ! grep -E "DataFrame|semType|ACL|viewer factory" src/panels/published-analysis-panel.ts | grep -v "^//" | grep -v "^import" | grep -v "DG\\.DataFrame" | grep -v "DG\\.Widget"; }</automated>
  </verify>
  <done>Panel exports `publishedAnalysisPanel(proteinId)`; first-line `isPublished(df)` guard; reads metadata column-first via `getPublishedMetadata`; renders 7 audit fields (Target, Shared, Shared by, Method, Comparison, FC cutoff, p cutoff); imports `buildMailtoUrl` from `../publishing/publish-state` (B-1 — NOT from share-dialog); user-visible strings pass jargon audit; strict TypeScript compiles.</done>
</task>

<task type="auto" priority="P2" tdd="false">
  <name>Task 2 (P2 - PUB-13 mailto sub-feature): Add supersede link + mailto button + edge cases</name>
  <files>src/panels/published-analysis-panel.ts</files>
  <read_first>
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md (D-04 "Newer version available" + PUB-13 mailto)
    - @packages/Proteomics/src/publishing/publish-state.ts (buildMailtoUrl signature — Plan 01; PUBLISHED_TAGS.SUPERSEDED_BY, PUBLISHED_TAGS.PUBLISHED_BY_EMAIL, META_COLUMNS parallels)
  </read_first>
  <action>
**B-3 task-level priority:** This task adds two things: the supersede link (P1 — PUB-07 audit context) AND the mailto button (P2 — PUB-13). The supersede link is non-deferrable; the mailto button can ship later. They share the same task because they both run after the 7 audit fields.

Add to `publishedAnalysisPanel` (after the 7 audit fields):

**Supersede link (D-04, conditional — P1):**
- Read `meta.supersededBy` from `getPublishedMetadata` (column FIRST per Pitfall 3 — already handled by Plan 01 helper). If non-null:
  - Build the link via `ui.link('Newer version available — open it', async () => { const newer = await grok.dapi.projects.find(meta.supersededBy!); await newer.open(); });`
  - Add the link at the TOP of the body (before the 7 audit fields) so the reviewer notices.
  - Per CONTEXT.md D-04: biologist's stale bookmark of v1 surfaces this link instead of forcing them to find v2 manually.
  - W-8 note: `getPublishedMetadata` reads `_meta_superseded_by` column FIRST, then `proteomics.superseded_by` tag. Plan 04 step 9 writes BOTH (W-8 dual-write); this panel automatically gets either path.

**Mailto button (PUB-13, P2 sub-feature):**
- Compute project name. Best source: `df.name` (Plan 02 set it to `<source>_published_<date>`) or `grok.shell.tv.parentView?.name` or via project lookup. For simplicity + reliability: use `df.name`.
- Read sharer email: column-first via `df.col(META_COLUMNS.PUBLISHED_BY_EMAIL)?.get(0)` then `df.getTag(PUBLISHED_TAGS.PUBLISHED_BY_EMAIL)`. May be null per RESEARCH Open Question 5.
- `buildMailtoUrl` is ALREADY IMPORTED in Task 1 from `../publishing/publish-state` (B-1 fix). Do NOT re-import from `../publishing/share-dialog`.
- Build URL:
  ```
  const mailtoUrl = buildMailtoUrl({
    sharerEmail: sharerEmail,
    sharerName: meta.publishedBy,
    projectName: df.name,
    publishedDateStr: meta.publishedAt instanceof Date ? meta.publishedAt.toISOString().slice(0,10) : String(meta.publishedAt).slice(0,10),
  });
  ```
- Render as `ui.link('Request re-run with different parameters', mailtoUrl)` — anchor with href, not callback (PUB-13: plain mailto: URL).
- Append at the BOTTOM of the body.

**Defensive edge cases:**
- If `df.col(META_COLUMNS.PUBLISHED_TARGET)` is null AND `df.getTag(PUBLISHED_TAGS.PUBLISHED_TARGET)` is null: the artifact is corrupt — render only `ui.divText('Shared analysis metadata unavailable. Ask the sharer to re-share.')`. This catches the case where both belt-and-braces sources failed (extremely defensive — should never happen if Plan 02 + Plan 04 ran correctly).
- If `meta.publishedAt` parse fails (corrupt date): render `'Shared: unknown date'` rather than crash.

JSDoc the function with: "Reviewer-side audit context panel. Reads metadata column FIRST, tag SECOND per Pitfall 3. Returns empty widget on non-published DataFrames. All user-facing strings pass Pitfall 14 jargon audit. Imports `buildMailtoUrl` from publish-state.ts (B-1 — avoids wave 2 -> wave 4 dependency cycle)."
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/panels/published-analysis-panel.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "supersededBy" src/panels/published-analysis-panel.ts | grep -qv '^0$' &amp;&amp; grep -c "buildMailtoUrl" src/panels/published-analysis-panel.ts | grep -qv '^0$' &amp;&amp; grep -c "Request re-run" src/panels/published-analysis-panel.ts | grep -qv '^0$' &amp;&amp; grep -c "Newer version available" src/panels/published-analysis-panel.ts | grep -qv '^0$' &amp;&amp; { ! grep "from '../publishing/share-dialog'" src/panels/published-analysis-panel.ts; }</automated>
  </verify>
  <done>Panel renders supersede link at top when `meta.supersededBy` present (P1 — uses W-8 dual-source read); renders mailto button at bottom using `buildMailtoUrl` imported from `../publishing/publish-state` (P2 — B-1 fix); corrupt-artifact fallback message; strict TypeScript compiles; no import from share-dialog.ts (B-1 verified).</done>
</task>

</tasks>

<threat_model>
## Trust Boundaries

| Boundary | Description |
|----------|-------------|
| reopened published DataFrame -> panel renderer | DataFrame may have stripped tags (Pitfall 3); panel must tolerate missing fields gracefully |
| reviewer's click on supersede link -> grok.dapi.projects.find | Server-side ACL enforces whether reviewer can see the newer version; panel surfaces error if find fails |

## STRIDE Threat Register

| Threat ID | Category | Component | Disposition | Mitigation Plan |
|-----------|----------|-----------|-------------|-----------------|
| T-15-02 | Information disclosure (panel renders on non-published live DF) | publishedAnalysisPanel | mitigate | First-line `isPublished(df)` guard returns empty widget — non-published DataFrames never see panel content |
| T-15-SP-08 | Tampering (mailto-link injection via metadata) | sharerEmail / meta.publishedBy → buildMailtoUrl (in publish-state.ts) | mitigate | `buildMailtoUrl` (Plan 01) wraps every variable in `encodeURIComponent`; metadata values flow through that encoder |
| T-15-SP-09 | Denial of service (panel crashes on corrupt metadata) | publishedAnalysisPanel | mitigate | Fallback messages for missing target/date; never throws, always returns a widget |
</threat_model>

<verification>
- TypeScript strict-mode compiles
- `tsc --noEmit src/panels/published-analysis-panel.ts` resolves with `buildMailtoUrl` imported from `../publishing/publish-state` (B-1)
- Panel returns empty widget on non-published DF (isPublished guard)
- Panel reads `META_COLUMNS.*` before `PUBLISHED_TAGS.*` (Pitfall 3 read order)
- Supersede link appears only when `supersededBy` is non-null
- Mailto link present when sharer email is non-null; absent (or no-recipient mailto) when null
- All user-facing strings pass Pitfall 14 jargon audit (manual review at code-review time)
- Plan-level priority is NOT P2 (B-3 — plan covers P1 requirements; only Task 2 mailto sub-feature for PUB-13 is P2)
</verification>

<success_criteria>
- `src/panels/published-analysis-panel.ts` exports `publishedAnalysisPanel(proteinId): DG.Widget`
- Panel renders 7 required audit fields (PUB-07)
- Supersede link wired to `grok.dapi.projects.find(supersededBy).open()` (P1)
- Mailto button uses Plan 01's `buildMailtoUrl` helper imported from `../publishing/publish-state` (P2; B-1)
- First-line `isPublished(df)` guard non-negotiable
- Corrupt-artifact fallback never crashes
- No import from `../publishing/share-dialog` (B-1 verified — would cause wave-2 → wave-4 compile failure)
</success_criteria>

<output>
Create `.planning/phases/15-read-only-publishing-foundation/15-06-SUMMARY.md` when done with:
- All 7 audit fields rendered, mapped to PUB-07
- Confirmation of column-first read order (Pitfall 3 mitigation)
- Notes on biologist-jargon-audit results (any strings flagged + replaced)
- Confirmation that `buildMailtoUrl` is imported from `../publishing/publish-state` (B-1)
- Confirmation that plan-level priority is removed; only mailto sub-task is P2 (B-3)
</output>

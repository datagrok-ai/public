# Polish plan — skills & artefacts to ship

The goal is **reliability**, not capability. Grokky's `Execute` provider can already do most of what the scenarios need by having the LLM write JS. The problem is non-determinism: the LLM has to remember API surface mid-call. Adding structured skills converts brittle generation into deterministic invocation.

This plan is organized as **classes of action**, not per-scenario cookie cutters. Each Tier 1 skill handles a whole category of audience prompts — "change the axis font," "different calculated column," "show me activity cliffs" — not just the exact prompts in the recorded demos.

## Two execution paths the agent already has

- **Registered platform functions** (Chem, Bio, anything with `@grok.decorators.func`) → invoked via MCP `call_function`. Signatures are discoverable.
- **js-api methods on objects** (`viewer.setOptions`, `df.selection.set`, `tv.getFiltersGroup().updateOrAdd`, …) → invoked via `datagrok-exec` blocks (LLM writes JS).

The polish skills below are the typed layer on top of both. Where a wrapper exists, the LLM calls the skill by name with validated params; where it doesn't, the LLM is improvising — that's the failure surface.

---

## Tier 1 — Generic class-of-action skills

Ten skills cover virtually every demo prompt — recorded or live. Each handles its whole category.

### T1.1 — `addViewer` / `configureViewer` (any viewer, any property)

**Signatures:**
- `addViewer(type: string, options: object): Viewer`
- `configureViewer(viewer: Viewer, options: object): void`

**Backing APIs:** `grok.shell.addViewer(type, options)`; `viewer.setOptions(options)`.

**Why a wrapper:** the property schema is large and the LLM misspells common ones. The wrapper validates `options` keys against the viewer's known properties (`viewer.props` reflects them) and gives a useful error instead of silent no-op.

**Audit-confirmed gotchas to validate against:**
- `showRegressionLine` (NOT `regressionLine`) — `js-api/src/interfaces/d4.ts:579`.
- Viewer type strings are exact and case-sensitive: `'Scatter plot'`, `'Trellis plot'`, `'Chem Similarity Search'`, etc. — `js-api/src/const.ts:691`.

**Covers live Q&A like:** "change the X-axis font," "color by a different column," "make it logarithmic," "swap to a 3D scatter," "hide the legend," "set a different aspect ratio."

**Backs scenarios:** S1, S2, S3, S4.

---

### T1.2 — `addCalculatedColumn` (any formula)

**Signature:** `addCalculatedColumn(df, name, formula, type?='auto'): Promise<Column>`

**Backing API:** `df.columns.addNewCalculated(name, expression, type?, treatAsString?, subscribeOnChanges?)` — `js-api/src/dataframe/column-list.ts:180`.

**DSL:** `${ColumnName}` row references, operators (`+ - * / ^`) and function calls (`Sub`, `Div`, `Round`, `HeavyAtomCount`, …) per `help/transform/add-new-column.md`. Both forms work; wrapper accepts either.

**Recompute semantics (audit-corrected):** formula-metadata changes trigger automatic recomputation; cell-value edits do not always — `subscribeOnChanges` governs source-column dependency. Demo claims must describe formula recompute, not cell-edit recompute.

**Covers live Q&A like:** "compute heavy atom count," "add LipE," "ratio of two columns," "log of activity," "any expression the DSL parses."

**Backs scenarios:** S2, S4.

---

### T1.3 — Table operations (any column, any criterion)

**Signatures:**
- `filterRows(view, column, criteria: {min?, max?, values?, substructure?})`
- `sortRows(df, columns: string[], orders?: ('asc'|'desc')[])`
- `selectRows(df, predicate: BitSet | ((i) => boolean) | indices: number[])`
- `aggregateBy(df, groupCols, aggMap: {[col]: 'sum'|'avg'|'count'|...})`
- `pivot(df, {rows, cols, value, agg})` / `unpivot(df, {idCols, valueCols})`
- `joinTables(left, right, {leftKey, rightKey, joinType, normalizeKey?: 'canonicalSmiles' | 'lowercase' | ...})`

**Backing APIs:** `tv.getFiltersGroup().updateOrAdd({...})` for filters — same internal pattern `Chem:Substructure Search` uses at `packages/Chem/src/package.ts:727`. `df.selection` BitSet at `data-frame.ts:161`. DataFrame aggregation / pivot / join via standard js-api methods.

**Generalization note:** `joinTables` with `normalizeKey` covers **catalog dedup** (chemistry: `canonicalSmiles`), biologics registries (sequence hash), and any join keyed on a normalized identifier. There is no separate dedup skill. *Implementation:* the wrapper accepts `normalizeKey`, adds a canonical column to both sides via the appropriate normalizer (e.g. `Chem:canonicalize` for SMILES), and invokes the underlying DataFrame join on the canonical column. Not a change to the base join API — it's a small preprocessing wrapper.

**Covers live Q&A like:** "filter to non-toxic," "sort by potency," "select the top 20," "group by target," "pivot by activity class," "join with the safety panel," "dedupe against our catalog."

**Backs scenarios:** S2, S3, S4 (incl. dedup), S5 (post-ACL filtering).

---

### T1.4 — `callRegisteredFunction` (any platform function)

**Signature:** `callRegisteredFunction(name: string, params: object): Promise<any>`

**Backing APIs:** `grok.functions.find({name})` + `func.apply(params)`, or MCP `call_function`. Function names and signatures are discoverable at runtime via `grok.functions.find` / MCP `list_functions`.

**Why this is one skill, not eleven:** every Chem function (`Chem:Substructure Search`, `Chem:Chem Space`, `Chem:R-Groups Analysis`, `Chem:Activity Cliffs`, `Chem:Scaffold Tree`, `Chem:Matched Molecular Pairs`, `Chem:Toxicity Risks`, …), every Bio function, every package's `@grok.decorators.func` registration is invoked the same way. The LLM doesn't memorize names — it discovers them, validates param shape against the registered signature, and invokes. New packages installed on a customer server become callable without changing Grokky.

**Verified-name reference table** lives in the appendix below — that's documentation for the audit trail, not eleven skills to build.

**Covers live Q&A like:** "show activity cliffs," "cluster by Morgan FP," "do a synthon search," "fit a curve," "anything any installed package registers."

**Backs scenarios:** S1, S3, S4, and any future scenario that touches a registered function.

---

### T1.5 — `createSharedProject` (project + views + share)

**Signature:** `createSharedProject({name, description?, viewsToSave?: SavedView[], shares?: {group, access: 'View'|'Edit'}[]}): Promise<{projectId, url}>`

**Backing APIs (audit-verified state):**
- ✅ MCP `create_project(name, description)` exists at `packages/Grokky/dockerfiles/mcp-server/src/index.ts:150`.
- ❌ No MCP `share_entity` tool today. Owned by **T1.6** (ships only in its gated, two-phase form). Until T1.6 lands, sharing falls back to `grok.dapi.permissions.grant(project, group, 'View')` via `datagrok-exec` — which T1.6 L2 will then block, so this fallback is transient by design.
- ❌ No MCP `save_layout` / `add_view_to_project` tool. Falls back similarly.
- ✅ `grok.dapi.permissions` and `grok.dapi.projects` exist at `js-api/src/dapi.ts:138`.

**Covers live Q&A like:** "save this as a dashboard for my team," "share with toxicology with view-only," "package these three views into a project."

**Backs scenarios:** S1, S4, S5.

---

### T1.6 — Mechanical confirmation gate for irreversible actions

**The gate is in code, not in instructions.** Three layers; each closes one of the loopholes that would let the LLM bypass confirmation.

#### L1 — Runtime-side interception of irreversible MCP tool calls

Define an `IRREVERSIBLE` tool set: `share_entity`, `delete_entity`, `overwrite_project`, future additions. **The MCP tool definitions are stubs** — text-returning placeholders, exactly like the existing `db_list_catalogs` / `db_describe_tables` / `db_try_sql` family at `packages/Grokky/dockerfiles/mcp-server/src/index.ts:281-331`. The actual logic lives in the claude-runtime, which intercepts these tool calls before they execute.

**Why stubs, not a server-side protocol:** the MCP contract returns `ToolResult = {content: [{type: 'text', text: string}]}` (`mcp-server/src/index.ts:10`). It has no native "needs confirmation" response shape. The `db_*` tools already solved this — they're registered for discoverability and intercepted in the runtime. Use the same pattern.

**The flow:**

1. Agent calls `share_entity({entityId, group, access})` via MCP.
2. Claude-runtime intercepts the call (same hook that handles `db_*`). The MCP stub never actually executes.
3. Runtime emits an `input_request` WebSocket event to the browser. Payload uses the existing `InputRequestEvent` shape (`packages/Grokky/src/claude/runtime-client.ts:17`) with `input: {kind: 'confirmation', action, summary, requestId}`.
4. Browser renders the confirmation card (L3). User clicks **Confirm** or **Cancel**.
5. Browser sends `respondToInput(sessionId, {confirmed: true | false, requestId})` via the existing channel (`runtime-client.ts:167-171`).
6. Runtime validates the `requestId` was one *it* issued for this session, and only then performs the actual share via `grok.dapi.permissions.grant(...)` (browser-side, also over the WebSocket).
7. Runtime returns the result as the MCP tool's return content.

**Why this is a gate:** the runtime is the only entity that knows whether the user confirmed. The LLM can't fabricate a `{confirmed: true}` response — the runtime ignores anything not flowing through `respondToInput` for a matching `requestId` it just issued. The MCP stub has no execution path at all; without runtime interception it returns the placeholder text and does nothing.

**State lives in the runtime, not the MCP server.** `mcp-server/src/index.ts:340-350` creates a fresh `McpServer` per HTTP request — no place to keep session state. The runtime container holds the pending-confirmation map keyed on `(sessionId, requestId)`.

#### L2 — `datagrok-exec` sandbox (close the side door)

Without this, the agent could write `grok.dapi.permissions.grant(...)` directly inside a `datagrok-exec` block and bypass L1 entirely.

**Today** (`packages/Grokky/src/claude/exec-blocks.ts:13-15`), the real SDK modules are passed straight into a `Function` constructor:

```ts
new Function('grok', 'ui', 'DG', 'view', 't',
  'return (async () => {' + code + '})()',
)(grok, ui, DG, view, t);
```

**Fix:** wrap `grok` (and `ui` if any irreversible methods live there) in a Proxy before passing in. The Proxy intercepts property paths that resolve to irreversible methods (`dapi.permissions.grant`, `dapi.projects.delete`, `dapi.scripts.delete`, etc.) and either:
- Throws with a redirect message ("use the `share_entity` MCP tool"), forcing the agent to take the gated path; or
- Routes the call through the same L1 `input_request` confirmation flow, so direct-dapi usage also lands at the gate rather than being blocked.

Either choice leaves the agent with **no callable** for the irreversible action outside the L1-gated path.

#### L3 — Confirmation message style (visible, not informational)

The `InputRequestEvent` channel is already plumbed end-to-end (`runtime-client.ts:17` defines the event; `:108-110` dispatches it to the `onInputRequest` Subject; `:167-171` is the response path). The work is purely on the browser-side renderer.

Extend the panel's `onInputRequest` handler (in `packages/Grokky/src/ai/panel.ts` or wherever the subscription lives) to recognize `input.kind === 'confirmation'` and render a distinct **action confirmation** card: border, icon, action summary in bold, **Confirm** / **Cancel** buttons. Buttons emit `respondToInput(sessionId, {confirmed: true | false, requestId})`.

This also delivers the inline-buttons UX we earlier deferred for the safety-net pattern — it lives here naturally.

#### Why this is a real gate

- **L1:** server refuses without a user-issued token.
- **L2:** no alternative API path the agent can reach.
- **L3:** user notices and acts.

The LLM cannot drift, skip, or be talked out of it — the code refuses. Pair with a soft system-prompt instruction as defense-in-depth so the agent volunteers the question rather than hits the wall, but **the wall is what makes "gate" honest in talking points.**

#### Earlier-rejected approaches

Recorded so they don't come back: structured `clarify()` skill with per-option metadata; trigger taxonomy; `package-knowledge.yaml` enrichment with `typicalQuestions` / `bestFor`; **system-prompt-only "soft form" as the primary mechanism** — insufficient because it's a convention, not a gate, and S5's IT-lead talking point becomes dishonest under it.

#### Definition of done

- The share MCP tool refuses to fire without a valid token, demonstrated by an attempted bypass that fails with `requires_confirmation`.
- `datagrok-exec` blocks attempting `grok.dapi.permissions.grant` throw or redirect.
- S5 records cleanly with the API-level refusal as the load-bearing beat.

#### Dependencies

`respondToInput()` / `input_request` channel (exists). MCP tool refactor (new code). `exec-blocks.ts` globals refactor (new code). Chat panel render addition (new code).

#### Backs scenarios

S5 (now honestly a gate). Bolsters S1 and S4 sharing steps — same mechanism applies wherever the agent reaches an irreversible MCP tool.

---

## Tier 2 — Convenience wrappers (only where the registered function has parameter friction)

These are not class-of-action skills. They exist because a specific registered function takes the *wrong* shape of input for what the LLM has.

### T2.1 — `injectSubstructureFilter(molColumn, smarts)`

**Why:** `Chem:Substructure Search` (`packages/Chem/src/package.ts:724`) opens an empty filter dialog with `molBlock: DG.WHITE_MOLBLOCK`. The LLM has a SMARTS string, not a UI to drive. This wrapper does what the function does internally (line 727 — `tv.getFiltersGroup().updateOrAdd({type: SUBSTRUCTURE, column, molBlock})`) but accepts SMARTS and parses it.

**Backs:** S3.

### T2.2 — `runRGroupsAnalysisHeadless(molColumn, {core: smarts | 'MCS'})`

**Why:** `Chem:R-Groups Analysis` (`:1095`) opens a dialog and waits for a user to click MCS + OK. Headless wrapper accepts the core (or `'MCS'` for auto-detect) and produces the R-columns + Trellis Plot in one call.

**Backs:** S3 (cleanest path), S4 (if R-groups added to the playbook).

### T2.3 — `setTrellisInnerViewer(viewer, innerType, valueColumn?)`

**Why (provisional):** the audit found trellis inner-viewer switching is event-driven (`d4-trellis-plot-viewer-type-changed`) AND requires Context Panel manipulation to set Value. If `viewer.setOptions({inner: {type, value}})` works as a single call, this collapses into T1.1 and we don't need it. **Verify before building.**

**Backs:** S3 step 4.

That's it for Tier 2. Anything not on this list goes through T1.4.

---

## Tier 3 — Platform / UX changes (not skills — backend or UI work)

### T3.1 — Partial-result signal in query results

**What:** Server-side change to surface `filtered_by_acl: true` (or pre-filter row count) in query result metadata.

**Why:** S5's "Bob got 1 of 4 trials" message requires the agent to know rows were filtered. Audit confirmed no such signal exists today.

**Backs:** S5 (the explicit graceful-degradation message; without this, the presenter voices it).

### T3.2 — Audit-log surfacing in chat

**What:** A clickable "audit" badge on agent messages, linking to the audit-log entry.

**Why:** Audit log already captures everything (`help/govern/audit/audit.md:148`); chat-side affordance is missing. UI work in `packages/Grokky/src/ai/`.

**Backs:** S5 (the visible-governance moment).

### T3.3 — Drag-drop chat-panel handler

**What:** DOM drop listener on the chat panel; uploads to `agents/` and triggers the existing file sync.

**Why:** File sync into `agents/` is verified working (`packages/Grokky/src/package.ts:62-75`); no drop handler on the chat itself. Today users drag into the File Manager.

**Backs:** S4 (the stretch beat; paste-the-md is the fallback).

### T3.4 — Remaining MCP tools

**What:** Add to `packages/Grokky/dockerfiles/mcp-server/src/`:
- `save_project_layout(projectId, layoutName, layoutJson)` (or `add_view_to_project(...)`)

**Note on `share_entity`:** previously listed here, now owned by **T1.6** as a gated tool (two-phase confirmation protocol). The bare `share_entity` is not added; it ships only in its gated form.

**Why:** T1.5 (`createSharedProject`) currently falls back to a `datagrok-exec` block calling `grok.dapi.projects.save` because the save-layout tool doesn't exist. With it, T1.5 is a clean single-shot composing `create_project` + `save_project_layout` + `share_entity` (T1.6).

**Backs:** S1, S4, S5 (cleans up the save-view path; sharing is T1.6).

---

## Tier 4 — Demo prep (one-time, not platform features)

### T4.1 — `grok s` demo-account preset script

**What:** Idempotent bash/`grok s` script that stages `alice@demo`, `bob@demo`, and the groups S5 expects (`clinical-research-all`, `project-alpha-only`, `clinical-ops`, `toxicology-review`) on a fresh demo server.

**Verified commands** per `tools/GROK_S.md:18-22`: `grok s users save`, `grok s groups save`, `grok s groups add-members`.

**Backs:** S5, S4 (`toxicology-review`).

### T4.2 — Curated Chembl query with `pChEMBLMin` parameter

**What:** Add to `packages/Chembl/queries/` — a parameterized query that filters bioactivity by `pchembl_value >= @pChEMBLMin`. Either a sibling of `activityDetailsForTarget.sql` or an extension.

**Why:** Audit confirmed `activityDetailsForTarget.sql` has no pChEMBL column or filter. S1's demo prompt requires the threshold.

**Backs:** S1.

### T4.3 — Pre-warm + pre-record harness

**What:** Script that runs each scenario end-to-end on the demo server, warms caches, validates each step, and captures clean fallback recordings.

**Backs:** every scenario; the final gate before the conference.

---

## Scenario coverage map

What each scenario actually depends on, in the new tier framing:

| Scenario | Tier 1 skills | Tier 2 wrappers | Tier 3 platform | Tier 4 demo prep |
|---|---|---|---|---|
| **S1** English → dashboard | T1.4 (`Chem:Chem Space`, Chembl query), T1.5 (project+share) | — | T3.4 (sharing MCP cleanup) | T4.2 (pChEMBL query) |
| **S2** Tabular transform | T1.2 (calc col), T1.3 (filter, select), T1.1 (scatter config) | — | — | — |
| **S3** Substructure SAR | T1.4 (Chem registered fns) | T2.1, T2.2, T2.3 | — | — |
| **S4** Drop the playbook | T1.3 (join for dedup), T1.4 (Chem registered fns), T1.5 (project) | T2.2 (if R-groups added) | T3.3 (drag-drop), T3.4 (sharing MCP) | — |
| **S5** Governance | T1.4 (curated queries), T1.5 (project+share), T1.6 (safety-net) | — | T3.1 (partial-result signal), T3.2 (audit badge) | T4.1 (demo accounts) |

---

## Scope recommendation

**Must-keep (irreducible reliability layer):** T1.1, T1.2, T1.3, T1.4, T1.5, T1.6. Six skills. Without these, scenarios become voiceover-with-handwaving.

**Strongly want:** T3.1 (S5's explicit beat depends on it), T4.1 (S5 needs the demo accounts), T4.2 (S1 needs the threshold query).

**Nice to have:** T2.x (fall back to T1.4 + manual UI clicks), T3.2 (presenter can voice the audit point), T3.3 (paste-the-md fallback), T3.4 (T1.5 works without it, just uglier), T4.3 (manual takes work).

**Cut order if scope slips** — least painful first:
T3.2 → T4.3 → T2.x → T3.3 → T3.4. **Never cut the Tier 1 six.**

T1.6 stays in must-keep despite being tiny — one-paragraph system-prompt change, zero new code paths.

---

## Live-demo Q&A coverage test

If the audience asks something off-script:

| Audience asks | Handled by |
|---|---|
| "Change the X-axis font" | T1.1 — `configureViewer(viewer, {xAxisFont: 'Arial 14'})` |
| "Different formula — heavy atom count" | T1.2 — `addCalculatedColumn(df, 'HAC', 'HeavyAtomCount(${molecule})')` |
| "Show me activity cliffs" | T1.4 — `callRegisteredFunction('Chem:Activity Cliffs', {...})` |
| "Filter to non-toxic" | T1.3 — `filterRows(view, 'ToxRisk', {values: ['none', 'low']})` |
| "Cluster by Morgan FP" | T1.4 — `callRegisteredFunction('Chem:BitBIRCH Clustering', {...})` |
| "Join with the safety panel" | T1.3 — `joinTables(...)` |
| "Dedupe against our catalog" | T1.3 — `joinTables(..., {normalizeKey: 'canonicalSmiles'})` |
| "Build a 3D scatter instead" | T1.1 — `addViewer('3D scatter plot', {...})` |
| "Save this for my team" | T1.5 — `createSharedProject({name, shares})` |
| "Sequence-based search instead" | T1.4 — `callRegisteredFunction('Bio:…', {...})` if Bio is installed |

Everything routes through ~10 generic skills. The scenarios just exercise specific combinations.

---

## Appendix — Verified registered names (documentation, not a build list)

This table is here so anyone reading a scenario can verify a cited name. T1.4 (`callRegisteredFunction`) handles every entry the same way; no per-row skill needed.

| Verb | Registered name | Implementation | File:line |
|---|---|---|---|
| Substructure search | `Chem:Substructure Search` | `SubstructureSearchTopMenu` | `packages/Chem/src/package.ts:724` |
| Similarity search | `Chem:Similarity Search` | `similaritySearchTopMenu` | `:548` |
| Diversity search | `Chem:Diversity Search` | `diversitySearchViewer` | `:557` |
| R-Groups Analysis | `Chem:R-Groups Analysis` | `rGroupsAnalysisMenu` | `:1095` |
| Chemical Space | `Chem:Chem Space` | `chemSpaceTopMenu` | `:878` |
| Canonicalize | `Chem:canonicalize` | | `:2513` |
| Strip salts/water | `Chem:removeWaterAndSalts` | | `:3010` |
| Structural Alerts | `Chem:Structural Alerts` | `structuralAlertsTopMenu` | top-menu `:1340` |
| Toxicity Risks | `Chem:Toxicity Risks` | | `widgets/toxicity.ts:58` |
| Drug Likeness | `Chem:Drug Likeness` | | `widgets/drug-likeness.ts:17` |
| Activity Cliffs | `Chem:Activity Cliffs` | | top-menu `:1156` |
| Matched Molecular Pairs | `Chem:Matched Molecular Pairs` | | top-menu `:2329` |
| Scaffold Tree | `Chem:Scaffold Tree` | | top-menu `:2291` |
| BitBIRCH Clustering | `Chem:BitBIRCH Clustering` | | top-menu `:748` |
| Chemical Properties | `Chem:Chemical Properties` | | top-menu `:2199` |

**Rule sets for `Chem:Structural Alerts`** (`packages/Chem/src/widgets/structural-alerts.ts:6-8`): `PAINS`, `BMS`, `SureChEMBL`, `MLSMR`, `Dundee`, `Inpharmatica`, `LINT`, `Glaxo`. Brenk is not included; add SMARTS to `files/alert-collection.csv` if needed.

# Phase 15: Read-Only Publishing Foundation — Research

**Researched:** 2026-06-07
**Domain:** `DG.Project` + `grok.dapi.permissions` + nested `grok.dapi.spaces`; deep-clone + tag-and-column belt-and-braces; round-trip verification + verify-and-rollback
**Confidence:** HIGH for stack and platform surface (everything verified in js-api at repo HEAD `21d37776ee` and against live server `release/1.27.3`); HIGH for D-03 (nested Spaces confirmed first-class in js-api AND in live server data); MEDIUM-HIGH for Pitfall 3 tag-survival shape (Phase 13 `e527d07ba1` direct evidence + js-api code surface, but a final enumeration script still needs to run at execution time on a real round-trip — Wave 0 builds it).

## Summary

Phase 15 is composition work on platform-stable primitives — there is no new math, no new state store, no new dependency, no schema. The shape is: trim DataFrame (deep-clone, allowlist 7 columns, write belt-and-braces metadata column) → save Project (single `addChild` per DataFrame, two DataFrames when enrichment is present) → grant View on the per-target child Space → verify-and-rollback if Edit slipped in → assert round-trip shape before returning success → optionally write `superseded_by` pointer on the prior version. Everything else (the dialog, the audit panel, the mailto, the slug, the version detection) is plumbing around those five steps.

The single load-bearing technical unknown — *do nested Datagrok Spaces exist as a first-class primitive* — is **resolved**: the live server returns Projects with namespace paths like `DDI:Target1:` (nested subspaces), and `js-api/src/dapi.ts:849` exposes `SpaceClient.addSubspace(child)` with explicit hierarchy support. The D-03 primary shape (umbrella Space contains per-target child Spaces) is implementable; the fallback shape is not needed.

**Primary recommendation:** Build `src/publishing/` exactly as laid out in CONTEXT.md §"Integration Points". Use `DG.Project.create()` + `project.addChild(tableInfo)` per the canonical 3-line ApiSamples pattern, extended to multiple DataFrames using the Bio test idiom (`project.addChild(tableInfo); project.addChild(layoutInfo)`). Grant View at the per-target Space level via `permissions.grant(space, group, false)`. Verify with `permissions.get(space)` (returns `{view: Group[], edit: Group[]}` — assert reviewer group is in `view` and NOT in `edit`). Round-trip assertion runs INSIDE `publishAnalysis` before returning success; failure rolls back via `dapi.projects.delete(project)`.

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Deep-clone DataFrame + write 7-column allowlist | Browser / Client (TypeScript) | — | Pure transformation on an in-memory `DG.DataFrame`; no server roundtrip until save |
| Write belt-and-braces metadata column | Browser / Client | — | Single-row column added to the cloned DF before `uploadDataFrame` |
| Persist DataFrame → Project | API / Datagrok server | Browser (orchestrates) | `grok.dapi.tables.uploadDataFrame` + `tables.save(tableInfo)` + `projects.save(project)` — server is the durable store |
| Ensure umbrella Space + per-target child Space exist | API / Datagrok server | Browser (idempotent check via `rootSpaceExists` / `subspaceExists`) | Spaces live on the server; the helper coordinates server calls |
| Grant View permission on Space | API / Datagrok server | Browser | `permissions.grant(space, group, false)` is a server-side ACL write |
| Verify-and-rollback gate | Browser (orchestrates) | API / Datagrok server (`permissions.get` + `projects.delete`) | Decision logic lives client-side; effects are server writes |
| Round-trip assertion (`assertPublishedShape`) | Browser | API / Datagrok server (`projects.find(id).open()`) | Client reads the reopened DF and asserts the contract |
| Audit panel (reviewer-side context) | Browser / Client | — | Pure rendering of tags + metadata column, no fetch |
| Mailto link | Browser / Client | — | `window.location.href = 'mailto:...'` or `<a href="mailto:...">` |
| Detect republish + write `superseded_by` | Browser | API / Datagrok server (`projects.filter` + `setTag` on prior Project) | Client constructs the filter; server runs it and persists the new tag |

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `datagrok-api` | `^1.25.0` (already installed; live server runs `1.27.3`) | `DG.Project`, `DG.DataFrame.clone`, `grok.dapi.{projects,tables,permissions,groups,spaces,users}`, `grok.user.current()`, `ui.dialog`, `ui.choiceInput`, `@grok.decorators.func`/`panel` | Platform-native, zero dependency cost, every primitive already used elsewhere in this codebase or in ApiSamples [VERIFIED: js-api/src/dapi.ts and js-api/src/entities/project.ts at HEAD] |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `@datagrok-libraries/test` | `^1.1.0` (already installed) | `category`, `test`, `expect`, `awaitCheck`, `delay` for `src/tests/publish-roundtrip.ts` | Round-trip test imports — same shape as `src/tests/volcano.ts` and Bio's `projects-tests.ts:14` |
| `cash-dom` | `^8.1.5` (already installed) | DOM manipulation inside audit panel | Already used in panels; do not reach for jQuery |
| `rxjs` | `^6.5.5` (already installed) | Subscription for cross-DF protein-highlight in published view (D-05) | Same pattern as v1.2 `enrichment-viewers.ts:9` |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Nested `DG.Space` (D-03 primary) | Project-nesting under one flat Space (D-03 fallback) | The fallback is unnecessary — nested Spaces are verified first-class (see §"Nested Spaces Verification" below). The fallback shape is documented in CONTEXT.md only because the researcher hadn't yet verified at discuss-phase time. |
| Multiple `addChild` calls (one per DF) | Single `addChild(layoutInfo)` only | Bio's `projects-tests.ts:37-38` adds BOTH `tableInfo` AND `layoutInfo` as children. For D-05's two DataFrames + one volcano layout, plan: `addChild(proteinTableInfo) + addChild(enrichmentTableInfo if present) + addChild(layoutInfo)`. Single-call alternatives drop either the data or the layout. |
| Per-target Space inheritance (D-03) | Per-Project grant only (skip Spaces) | Per-Project grant works (`permissions.grant(project, group, false)` per `dapi/layouts-and-permissions.js`), but loses the audit-by-target organization PUB-04 wants. Both versions of the same target end up loose in the user's namespace with no parent. Spaces give the per-target audit folder for free. |
| Hand-rolled `mailto:` builder | `ui.element('a', {href: 'mailto:...'})` direct | Mailto URLs are 100% client-side; no library needed. `encodeURIComponent` on subject/body is the only escaping required. |

**Installation:** No new dependencies. v1.3 `package.json` is unchanged for v1.4 (confirmed at `packages/Proteomics/package.json:14-25`).

**Version verification:**
- `datagrok-api ^1.25.0` — Confirmed at `packages/Proteomics/package.json:23`. Live server runs `release/1.27.3` per `grok s raw GET /api/info/server`. All Phase 15 APIs (`Project.addChild`, `permissions.grant`, `permissions.get`, `spaces.createRootSpace`, `SpaceClient.addSubspace`, `User.current()`, `dapi.projects.{save,find,delete,filter}`) are present at js-api HEAD `21d37776ee`. [VERIFIED: file reads]
- `@datagrok-libraries/test ^1.1.0` — already installed. [VERIFIED: package.json:22]

## Package Legitimacy Audit

**Not applicable for this phase.** Phase 15 installs ZERO new packages. Every primitive is already in `datagrok-api` (already-declared dep) or in existing `@datagrok-libraries/*` (already-declared deps). Slopcheck protocol is moot when the package list is empty.

| Package | Registry | Age | Downloads | Source Repo | slopcheck | Disposition |
|---------|----------|-----|-----------|-------------|-----------|-------------|
| *(none — no new packages)* | — | — | — | — | — | — |

## Architecture Patterns

### System Architecture Diagram

```
Proteomics | Share | Share Analysis for Review...           (NEW menu item)
  │
  ▼
showShareForReviewDialog(activeDf)                          (src/publishing/share-dialog.ts)
  ├── precondition gate: activeDf.getTag('proteomics.de_complete') === 'true'
  ├── findPriorShare(target, group) via projects.filter      (NEW: src/publishing/publish-state.ts)
  │    └── if found → pre-fill inputs + banner "Will publish as v<N+1>, supersede '<prior>'"
  ├── ui.dialog: target string + ChoiceInput<Group> + note textarea
  └── on OK ─────────► publishAnalysis(activeDf, {target, group, note, priorVersion?})
                          │
                          ▼
                  publishAnalysis()                          (src/publishing/publish-project.ts)
                          │
                          ├── 1. trimForPublish(activeDf)    (NEW: src/publishing/trim-dataframe.ts)
                          │     • df.clone(undefined, allowlist7)
                          │     • re-set every proteomics.* tag explicitly
                          │     • addMetadataColumn() — single-row typed columns:
                          │         published_target, published_at, published_by,
                          │         published_de_method, published_fc_threshold,
                          │         published_p_threshold, published_version, published_id
                          │     • frozen.name = `<source>_published_<date>`
                          │     • mark published / published_at / published_by tags
                          │
                          ├── 2. trimEnrichmentForPublish(activeDf?)  (opportunistic, D-05)
                          │     • if shell has enrich DF tagged proteomics.enrichment
                          │     • allowlist: Term Name, Source, FDR, p-value,
                          │       Intersection (+ Direction if Phase 13 split present)
                          │
                          ├── 3. ensureUmbrellaSpace('Proteomics-Reviews')   (D-03)
                          │     • spaces.rootSpaceExists('Proteomics-Reviews')
                          │     • if not: spaces.createRootSpace(...)
                          │
                          ├── 4. ensurePerTargetChildSpace(slug)               (D-03)
                          │     • umbrellaClient.subspaceExists(child)
                          │     • if not: umbrellaClient.addSubspace(child)
                          │
                          ├── 5. saveProject(trimmedProteinDf, enrichDf?)
                          │     • project = DG.Project.create()
                          │     • project.name = `Proteomics-Review-<slug>-v<N>-<date>`
                          │     • project.addChild(proteinTableInfo)
                          │     • project.addChild(enrichTableInfo) if present
                          │     • project.addChild(layoutInfo)
                          │     • await dapi.tables.uploadDataFrame(...) for each DF
                          │     • await dapi.tables.save(tableInfo) for each
                          │     • await dapi.views.save(layoutInfo)
                          │     • await dapi.projects.save(project)
                          │     • move project into child Space:
                          │       childSpaceClient.addEntity(project.id)
                          │
                          ├── 6. grantViewToReviewerGroup(childSpace, group)   (D-03)
                          │     • await dapi.permissions.grant(childSpace, group, false)
                          │
                          ├── 7. VERIFY-AND-ROLLBACK GATE                       (non-negotiable #1)
                          │     • perm = await dapi.permissions.get(childSpace)
                          │     • assert group is in perm.view AND NOT in perm.edit
                          │     • on failure: dapi.projects.delete(project); throw
                          │
                          ├── 8. assertPublishedShape ROUND-TRIP GATE           (non-negotiable #2)
                          │     • const fresh = await (await dapi.projects.find(project.id)).open()
                          │     • assert: every required proteomics.* tag present
                          │       OR readable from metadata column
                          │     • assert: PROTEIN_ID semType present
                          │     • assert: trimmed.name matches
                          │     • assert: volcano viewer round-tripped (Bio
                          │       projects-test pattern: search tv.viewers for type)
                          │     • on failure: dapi.projects.delete(project); throw
                          │
                          └── 9. UPDATE SUPERSEDE CHAIN if priorVersion        (D-04)
                                • priorProj.setTag('proteomics.superseded_by', project.id)
                                • newProj.setTag('proteomics.supersedes', priorProj.id)
                                • return project
```

### Recommended Project Structure

```
src/
├── publishing/                    # NEW directory (per ARCHITECTURE.md §"src/publishing/")
│   ├── publish-state.ts           # tag helpers + isPublished + findPriorShare + slug
│   ├── trim-dataframe.ts          # trimForPublish + trimEnrichmentForPublish
│   ├── publish-project.ts         # publishAnalysis (the orchestrator + both gates)
│   ├── share-dialog.ts            # showShareForReviewDialog
│   └── assert-published-shape.ts  # round-trip assertion helper
├── panels/
│   ├── uniprot-panel.ts           # unchanged
│   └── published-analysis-panel.ts # NEW: reviewer-side audit + "Request re-run" mailto
├── tests/
│   └── publish-roundtrip.ts       # NEW: load-bearing gate per success criterion 3
├── package.ts                     # +1 menu item, +1 panel decorator
└── utils/proteomics-types.ts      # NO new SEMTYPEs anticipated (direction may stay
                                    # string-typed without a SEMTYPE — Claude's discretion)
```

`detectors.js` — UNCHANGED for Phase 15 (no new SEMTYPEs).

### Pattern 1: Canonical `DG.Project` save + multi-DataFrame round-trip

**What:** The save-then-find-then-open cycle that the entire Phase 15 contract sits on.
**When to use:** Inside `publishAnalysis` step 5 AND inside `assertPublishedShape` step 8.
**Example:**
```typescript
// Source: packages/ApiSamples/scripts/dapi/projects.js (canonical 3-line shape)
//         + packages/Bio/src/tests/projects-tests.ts:26-47 (multi-child + layout + reopen)
const project = DG.Project.create();
project.name = `Proteomics-Review-${slug}-v${version}-${dateStr}`;

const proteinTableInfo = trimmedProteinDf.getTableInfo();
project.addChild(proteinTableInfo);
await grok.dapi.tables.uploadDataFrame(trimmedProteinDf);
await grok.dapi.tables.save(proteinTableInfo);

if (trimmedEnrichDf) {
  const enrichTableInfo = trimmedEnrichDf.getTableInfo();
  project.addChild(enrichTableInfo);
  await grok.dapi.tables.uploadDataFrame(trimmedEnrichDf);
  await grok.dapi.tables.save(enrichTableInfo);
}

const layoutInfo = volcanoTableView.getInfo();
project.addChild(layoutInfo);
await grok.dapi.views.save(layoutInfo);

await grok.dapi.projects.save(project);

// === round-trip verification (assertPublishedShape) ===
const projId = project.id;
const reopened = await grok.dapi.projects.find(projId);
await reopened.open();
// inspect: grok.shell.tv.dataFrame.getTag('proteomics.published') === 'true'
// inspect: grok.shell.tv.dataFrame.col('Protein ID')?.semType === SEMTYPE.PROTEIN_ID
// inspect: a volcano viewer is in tv.viewers
```

### Pattern 2: Nested-Space create + permission grant + verify

**What:** Umbrella Space → per-target child Space → reviewer-group View grant → readback.
**When to use:** Inside `publishAnalysis` steps 3-4 and step 7 (verify-and-rollback gate).
**Example:**
```typescript
// Source: packages/ApiSamples/scripts/dapi/spaces.js +
//         packages/ApiSamples/scripts/dapi/layouts-and-permissions.js +
//         js-api/src/dapi.ts:809-866

// Step 3 — ensure umbrella
let umbrella: DG.Project;  // Spaces ARE Projects with isSpace === true
const UMBRELLA_NAME = 'Proteomics-Reviews';
if (await grok.dapi.spaces.rootSpaceExists(UMBRELLA_NAME)) {
  umbrella = await grok.dapi.projects.filter(`name = "${UMBRELLA_NAME}" and isSpace = true`).first();
} else {
  umbrella = await grok.dapi.spaces.createRootSpace(UMBRELLA_NAME);
}
const umbrellaClient = grok.dapi.spaces.id(umbrella.id);

// Step 4 — ensure per-target child
const childName = `Proteomics-Review-${slug}`;
let childSpace: DG.Project;
if (await umbrellaClient.subspaceExists(childName)) {
  // No direct subspaceById helper — list children of umbrella, find by name
  const children = await umbrellaClient.children.filter('Project', false).list();
  childSpace = children.find((c) => c.friendlyName === childName) as DG.Project;
} else {
  childSpace = await umbrellaClient.addSubspace(childName);
}
const childClient = grok.dapi.spaces.id(childSpace.id);

// (Project save happens here, then moved into childSpace via childClient.addEntity(project.id))

// Step 6 — grant View to reviewer group at child-Space level
await grok.dapi.permissions.grant(childSpace, reviewerGroup, /* edit */ false);

// Step 7 — VERIFY-AND-ROLLBACK GATE
// permissions.get returns { view: Group[], edit: Group[] } per dapi.ts:682-687
const perm: {view: DG.Group[]; edit: DG.Group[]} = await grok.dapi.permissions.get(childSpace) as any;
const inView = perm.view.some((g) => g.id === reviewerGroup.id);
const inEdit = perm.edit.some((g) => g.id === reviewerGroup.id);
if (!inView || inEdit) {
  await grok.dapi.projects.delete(project);   // rollback the just-saved Project
  throw new Error(
    'Reviewer group already has elevated access via Space inheritance — publish aborted; ' +
    'ask an admin to scope the umbrella Space\'s permissions',
  );
}
```

### Pattern 3: Belt-and-braces metadata column

**What:** Encode every critical tag value as a single-row typed column on the trimmed DF.
**When to use:** Inside `trimForPublish`, after the clone and tag-set steps.
**Example:**
```typescript
// Source: Pitfall 3 mitigation; CONTEXT.md Claude's-discretion "Single-row metadata
//         column shape (PUB-11, Pitfall 3)"

function addMetadataColumns(frozen: DG.DataFrame, meta: PublishMetadata): void {
  // Single-row pattern: append one row, write each field into a 1-row typed column.
  // Reader (published-analysis-panel.ts) reads col(name).get(0).
  // Read column FIRST, tag SECOND (column survives serializer better per Phase 13 e527d07ba1).
  const r = frozen.rowCount;  // we add metadata *as columns*, not as a separate row
  frozen.columns.addNewString('_meta_published_target').init(() => meta.target);
  frozen.columns.addNewDateTime('_meta_published_at').init(() => meta.publishedAt);
  frozen.columns.addNewString('_meta_published_by').init(() => meta.publishedBy);
  frozen.columns.addNewString('_meta_published_de_method').init(() => meta.deMethod);
  frozen.columns.addNewFloat('_meta_published_fc_threshold').init(() => meta.fcThreshold);
  frozen.columns.addNewFloat('_meta_published_p_threshold').init(() => meta.pThreshold);
  frozen.columns.addNewInt('_meta_published_version').init(() => meta.version);
  frozen.columns.addNewString('_meta_published_id').init(() => meta.publishId);
  // Optional: hide these columns from default grid view via col.setTag('.color.coding', null)
  // and col.setTag('.hidden', 'true') so the reviewer doesn't see metadata as data.
}
```

**Note on column-vs-row encoding:** Each metadata field is its own column (every row gets the same value). This costs `N rows × 8 metadata fields × 8 bytes = trivial` (a 5,000-protein dataset adds ~320 KB). The alternative — appending a "metadata row" — breaks the row semantics ("each row is a protein"). Columns-with-constant-value preserve semantics and survive the serializer per Phase 13 evidence.

### Anti-Patterns to Avoid

- **`df.columns.remove(...)` on the live source DF** (Pitfall 1): the source DF is shared with the user's view. Mutating it during publish corrupts the analyst's working state AND leaks future re-runs into the "frozen" Project. ALWAYS clone first.
- **Skipping `assertPublishedShape` because "save succeeded"** (Pitfall 3): Phase 13 commit `e527d07ba1` proved that save-success does not guarantee reopen-fidelity. The round-trip gate runs INSIDE `publishAnalysis`, not just in CI.
- **Per-Project permission grant instead of per-Space** (D-03): per-Project works but doesn't audit-by-target (PUB-04 intent). Use per-Space grant so all versions for one target inherit.
- **Granting Edit via inheritance side-effect** (Pitfall 2): if the umbrella Space already grants Edit to a group, the child Space inherits. `permissions.get` MUST be checked AFTER grant and rollback issued if Edit slipped in.
- **Hard-coded column names** (CLAUDE.md convention): every column lookup goes through `findColumn(df, semType, nameHints)` or `findProteomicsColumns(df)` so the trim helper works against any vendor's parser output, not just MaxQuant's.
- **`fetch(url)` direct call** (CLAUDE.md, package convention): mailto: is the only external-link surface in Phase 15. No `fetch` calls. If a future requirement adds one, use `grok.dapi.fetchProxy`.
- **Storing per-field audit tags instead of a metadata column** (Anti-Pattern 2 in ARCHITECTURE.md): every metadata field also goes in the column for serializer survival.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Project serialization | Custom JSON of DataFrame + viewers + layout | `DG.Project.create() + addChild + dapi.tables.uploadDataFrame + dapi.projects.save` | Native; survives platform schema migrations; reopens via standard UI |
| ACL grant/check | Custom user/group membership map | `dapi.permissions.grant(entity, group, edit)` + `dapi.permissions.get(entity)` | Server enforces; client-side ACLs are advisory only |
| Project hierarchy | Folder-of-Projects naming scheme | `dapi.spaces.createRootSpace` + `SpaceClient.addSubspace` | Spaces give per-target organization for free + ACL inheritance (which we must verify, per D-03) |
| Group lookup / membership | Reading `groups.json` from AppData | `dapi.groups.list()` + `ChoiceInput` (D-02) | Real-time, server-authoritative, matches platform muscle memory |
| Slug sanitization | Multi-line regex chain | One regex + cap: `slug = target.replace(/[^A-Za-z0-9._-]+/g, '-').replace(/^-+|-+$/g, '').slice(0, 64)` | Stdlib regex is enough; "freeform → safe filename" is a tiny problem |
| Project deletion (rollback) | Manual cleanup of children | `dapi.projects.delete(project)` | Standard delete cascades to children per `HttpDataSource.delete` |
| Mailto link | Email service integration | Plain `mailto:` URL with `encodeURIComponent` on subject + body | PUB-13 is explicitly mailto, not SMTP |
| User friendly name | Custom user lookup | `grok.shell.user.friendlyName` (`Entity.friendlyName` per `entities/entity.ts:36`) — set server-side, not client-trustable | Never accept client-supplied `published_by` (security pitfall in PITFALLS.md "Security Mistakes" #7) |
| Version detection | Linear scan of project list | `dapi.projects.filter(\`name like "Proteomics-Review-${slug}-v%"\`).list()` then `max(version)+1` | Server-side smart filter — see §"Open Question 6 — Republish Detection Query" |
| Round-trip layout assertion | Reading viewer DOM | Iterate `tv.viewers` and assert by `.type` — Bio's `projects-tests.ts:65-73` pattern | Test framework convention; deterministic |

**Key insight:** Every piece of Phase 15 has a platform-native primitive. The work is composition, not invention. The two non-negotiable gates (`assertPublishedShape` and verify-and-rollback) are the entire risk surface; everything else is a 5-line call to `dapi.X.Y(...)`.

## Runtime State Inventory

> Phase 15 is a **greenfield additive** phase, not a rename/refactor. This section is normally omitted, but Phase 15 has one wrinkle: it INTERACTS with existing pipeline state on the source DataFrame (it reads tags, then clones). For completeness, here is the state Phase 15 reads, writes, or could regress:

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | None — Phase 15 introduces no new persisted store outside the standard `DG.Project` save path. AppData is NOT touched (per ARCHITECTURE.md §12 — Phase 17 is the first to touch AppData). | None |
| Live service config | None — no n8n, no external service config. The mailto link is a URL, not a service. | None |
| OS-registered state | None — no Task Scheduler, no pm2, no launchd touched. | None |
| Secrets / env vars | None — Phase 15 has zero credential handling (`packages/Proteomics/CLAUDE.md` line 13: "Auth / credentials: none"). | None |
| Build artifacts | None — no new SEMTYPEs (`direction` is staying as string per Claude's discretion in CONTEXT.md), so no `detectors.js` edit, so no stale-cached detector issue. | None |

**Source DataFrame mutation risk (Pitfall 1):** The single runtime-state risk is the source DF. Phase 15 reads tags + columns from it; if the implementation accidentally calls `df.columns.remove(...)` or `df.setTag(...)` on the live source instead of the clone, the analyst's working state is corrupted. Mitigation lives in the `trimForPublish` contract: every mutation is on the clone returned by `df.clone(undefined, allowlist)`. The Wave 0 test (`publish-roundtrip.ts`) explicitly asserts source-unchanged after publish.

## Common Pitfalls

### Pitfall 1: Stale-snapshot leak (Pitfall 1 in PITFALLS.md)
**What goes wrong:** Publish forgets to deep-clone. The "frozen" project shares the live DF reference; the analyst re-runs DE on the source later; the reviewer sees changed numbers.
**Why it happens:** Architectural default in this package is "mutate in place" (ARCHITECTURE.md §"In-place Mutation"); `df.columns.remove(...)` is a tempting one-liner that mutates not clones.
**How to avoid:** `publishAnalysis` first line: `const frozen = df.clone(undefined, allowlist7);`. Then `frozen.setTag(...)`, `frozen.columns.addNewString(...)`, etc. Never touch `df` directly.
**Warning signs:** Reviewer says "this looks different from last week"; row count > 7 columns after re-open; mutating one DF visibly affects the other in a debug session.

### Pitfall 2: Space-inheritance Edit leak (Pitfall 2 in PITFALLS.md; CONTEXT.md D-03 mitigation)
**What goes wrong:** Reviewer group has Edit on the umbrella Space (from prior unrelated grant). Child Space inherits. Reviewer can edit the published clone.
**Why it happens:** Datagrok's permission model is per-Entity AND per-Space; the per-Space grant doesn't shadow or revoke inherited permissions.
**How to avoid:** Step 7 of `publishAnalysis` is the verify-and-rollback gate (see Pattern 2 above). Read `permissions.get(childSpace)`; assert `view` contains reviewer group AND `edit` does NOT. If assertion fails, delete the new Project and surface the error.
**Warning signs:** Reviewer can right-click a column and see "Delete column"; the verify gate returns `inEdit === true`.

### Pitfall 3: Serializer strips tags / semTypes / df.name (Pitfall 3 in PITFALLS.md; commit `e527d07ba1`)
**What goes wrong:** Reviewer opens the project. Volcano is blank. `Display Name` column lost its `SEMTYPE.DISPLAY_NAME`. UniProt panel doesn't auto-fire. `df.name` is "table1".
**Why it happens:** Phase 13 commit `e527d07ba1` directly proves the platform serializer strips at least `look.filters[]` and `look.columnNames` regardless of input shape. The general contract is partial; the project saves successfully, but some tags/columns vanish on reopen.
**How to avoid:** TWO defenses, both required:
  1. **Belt-and-braces metadata column** — every critical tag value also lives in a single-row typed column (Pattern 3 above). Reader (panel) reads column FIRST, tag SECOND.
  2. **`assertPublishedShape` round-trip gate** — `publishAnalysis` step 8 actually performs the save → find → open cycle and asserts every required tag, semType, df.name, and viewer is present on the reopened DF. If anything missing, rollback.
**Warning signs:** Reopen test fails an assertion; reviewer sees blank volcano; `df.name` reads `table1`; `col.semType` is null on Protein ID.

### Pitfall 4: Versioning ambiguity (Pitfall 4 in PITFALLS.md; CONTEXT.md D-04)
**What goes wrong:** Republish overwrites or shadows a prior version. Reviewer's bookmark breaks; old numbers are gone; audit trail destroyed.
**How to avoid:** Republish creates a NEW Project with `+1` version number (`Proteomics-Review-<slug>-v<N+1>-<date>`). Soft pointer `proteomics.superseded_by = <new_id>` on the prior; bidirectional `proteomics.supersedes = <prior_id>` on the new. NEVER overwrite or delete the prior. The audit panel surfaces "Newer version available: [link]" when `superseded_by` is present.
**Warning signs:** Republish reuses the same project id; old URL returns 404; the project name contains a duplicate timestamp suggesting overwrite.

### Pitfall 14: Reviewer-side jargon / hung viewers (Pitfall 14 in PITFALLS.md; CONTEXT.md D-domain pin)
**What goes wrong:** Audience is biologist-consumer (CONTEXT.md domain pin). Words like "DataFrame", "tag", "semType", "ACL", "viewer factory" appear in the dialog or audit panel; biologist disengages.
**How to avoid:** Banned-word list applied to every reviewer-touchable string: dialog labels, panel text, mailto body, error messages. Replacements: "DataFrame" → "table"; "tag" → "label"; "ACL"/"permissions" → "who can see"; "viewer factory" → "view". Every async load on reviewer side has `DG.TaskBarProgressIndicator` with explicit phase messages (PITFALLS.md UX section item 2).
**Warning signs:** Dialog string contains a banned word; viewer renders >500ms with no spinner; published Project's title is just a UUID.

### Phase 15-specific pitfall: Volcano color lock + Filters viewer must survive into published view
**What goes wrong:** The Phase 14 D-04 magenta/cyan/gray volcano color lock and the Phase 14 D-05 Filters-viewer-scoped-to-Comparison-and-Gene-name-and-Protein-ID pattern were tuned for the analyst's view. If the trim drops the columns those features depend on, the published volcano will look different from the expert's working view (D-domain pin requires they match).
**Why it happens:** The 7-column allowlist (PUB-02) doesn't explicitly list Subcellular Location or Comparison columns. The volcano color binding silently falls back; the Filters viewer dock loses scope.
**How to avoid:** The trim allowlist for the published volcano must include any column the volcano color/filter bindings DEPEND ON. At minimum: Protein ID, Gene Name (Display Name), log2FC, p-value, adj.p-value, significant, direction. If `Subcellular Location` is bound (D-05 of Phase 13), include it too. If `Comparison` is bound for Spectronaut Candidates, include it. The audit panel's "Locked view" UX should match what the expert showed.
**Warning signs:** Reopened volcano shows gray dots only (direction column lost); volcano shows raw axis labels (Phase 14 G1 regression); Filters viewer is empty.

## Code Examples

Verified patterns from official sources:

### Save + open round-trip with layout (Bio pattern)

```typescript
// Source: packages/Bio/src/tests/projects-tests.ts:26-47
const project = DG.Project.create();
project.name = 'Test project';
const tableInfo = tv.dataFrame.getTableInfo();
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
```

### Grant View + read-back

```typescript
// Source: packages/ApiSamples/scripts/dapi/layouts-and-permissions.js
//         + js-api/src/dapi.ts:682-707
let adminUser = await grok.dapi.users.filter('login = "admin"').first();
await grok.dapi.permissions.grant(allLayouts[1], adminUser.group, true);  // true = edit
console.log(await grok.dapi.permissions.get(allLayouts[1]));
// Output shape (verified against dapi.ts:682): { view: [Group, ...], edit: [Group, ...] }
```

### Nested Space create + subspace + entity

```typescript
// Source: packages/ApiSamples/scripts/dapi/spaces.js:6-67
const exists = await grok.dapi.spaces.rootSpaceExists(spaceName);
const space = exists
  ? await grok.dapi.projects.filter(`name = "${spaceName}" and isSpace = true`).first()
  : await grok.dapi.spaces.createRootSpace(spaceName);
const client = grok.dapi.spaces.id(space.id);

// Nested subspace
const childSpace = await client.addSubspace('ChildSpace');

// Move a saved Project into the child Space
await grok.dapi.spaces.id(childSpace.id).addEntity(savedProject.id);
```

### Filter list of projects by tag (republish detection)

```typescript
// Source: js-api/src/dapi.ts:329 (filter — smart-filter expression)
//         + general HttpDataSource pattern
// Tag-based filtering: smart-filter supports `tag` predicates; exact syntax for
// custom (proteomics.*) tags is the open question — verify at plan time.
// Two candidate shapes:
//   const matches = await grok.dapi.projects.filter(
//     `name like "Proteomics-Review-${slug}-v%"`
//   ).list();
// — then client-side filter by reviewer-group access.
```

### `DG.DataFrame.clone` with column allowlist + tag reset

```typescript
// Source: src/viewers/qc-computations.ts pattern + heatmap.ts createExpressionHeatmap;
//         CLAUDE.md "Heatmap clones the DataFrame for filter isolation"
const allowlist = ['Protein ID', 'Gene Name', 'log2FC', 'p-value',
                   'adj.p-value', 'significant', 'direction'];
const frozen = df.clone(undefined, allowlist);
// IMPORTANT: clone tag-preservation is partial per Pitfall 3.
// Re-set every required tag explicitly:
frozen.setTag('proteomics.source', df.getTag('proteomics.source'));
frozen.setTag('proteomics.de_method', df.getTag('proteomics.de_method'));
frozen.setTag('proteomics.de_complete', 'true');
frozen.setTag('proteomics.groups', df.getTag('proteomics.groups'));
frozen.setTag('proteomics.published', 'true');
frozen.setTag('proteomics.published_at', isoDate);
frozen.setTag('proteomics.published_by', grok.shell.user.friendlyName);
frozen.setTag('proteomics.published_target', target);
frozen.setTag('proteomics.published_id', publishId);
frozen.setTag('proteomics.published_version', String(version));
// Then add belt-and-braces metadata columns (Pattern 3 above)
frozen.name = `${df.name}_published_${dateStr}`;
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Publish via Excel + email | `DG.Project` + permissions-grant + reviewer-opens-link | v1.4 (this phase) | First time the package supports cross-team review at all |
| Per-Project permission grant | Per-Space grant + inheritance | v1.4 (this phase D-03) | Audit-by-target organization; all versions of one target share an ACL surface |
| Tag-only metadata | Tag + single-row metadata column (belt-and-braces) | v1.4 (this phase Pitfall 3) | Survives serializer-strip; future-removable when platform fixes (per saved memory `feedback_keep_workaround_capture_future`) |
| Overwrite republish | Supersede chain (`superseded_by` + `supersedes`) | v1.4 (this phase D-04 / PUB-10) | Audit trail preserved; reviewer bookmarks stable |
| Trust client-supplied `published_by` | Read `grok.shell.user.friendlyName` server-authoritative | v1.4 (this phase, per PITFALLS.md Security Mistakes #7) | Audit log can't be forged |

**Deprecated / outdated:**
- **Project-nesting fallback (D-03 fallback)** — not needed. Nested Spaces are verified first-class in js-api `1.27.3` AND in live server data. Plan against the primary shape only.
- **Overwrite republish** — explicitly rejected per CONTEXT.md "Deferred Ideas" (loss of audit trail).

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | `grok.dapi.permissions.get(entity)` returns `{view: Group[], edit: Group[]}` exactly. | Pattern 2, Pitfall 2 mitigation | LOW — surface is documented at js-api/src/dapi.ts:682-687 explicitly; the `view` and `edit` field assignment is `toJs`-converted. But: the Map has additional fields (`'Share'`, `'Delete'`) per `permissions.check`'s typed parameter union. Verification at plan-build time: print `Object.keys(perm)` once and confirm. [VERIFIED: js-api/src/dapi.ts:682-687] |
| A2 | The reviewer-group's grant on the per-target child Space propagates to Projects subsequently moved INTO that Space (i.e., grant Space first, then `childClient.addEntity(project.id)` is sufficient — no per-Project regrant). | Pattern 2 step 6 + step 7 | MEDIUM — Datagrok docs say "When an entity is moved to a new project, it adopts the permissions of the new project" but no live API trace at research time. Wave 0 test (`publish-roundtrip.ts`) must call `permissions.get(reopenedProject)` AND `permissions.get(childSpace)` separately to confirm inheritance works as expected. If it doesn't, alternative is per-Project grant + skip Space grant. |
| A3 | `dapi.projects.delete(project)` cascades to the child `TableInfo` entities created by `addChild`. | Verify-and-rollback step 7 / round-trip step 8 | MEDIUM — `HttpDataSource.delete` is one-line `grok_DataSource_Delete` (dapi.ts:296). The cascade behavior depends on server. If it doesn't cascade, the rollback path may leave orphan `TableInfo` rows (cosmetic, not security). Wave 0 test should call `dapi.tables.find(tableInfoId)` after rollback and confirm 404 / NotFoundException. |
| A4 | Tags on a `DG.Project` (not on the DataFrame inside it) — specifically `proteomics.superseded_by` and `proteomics.supersedes` — survive `projects.save` and are readable via `project.getTag(...)` or `project.options['proteomics.superseded_by']`. | Pattern 1 + D-04 republish | MEDIUM — `Project.options` is a `MapProxy` (project.ts:26) — a key-value store on the Project, distinct from `DataFrame.tags`. Plan should use `project.options['proteomics.superseded_by'] = newId` rather than `project.setTag(...)` (which is a DataFrame method). Round-trip test asserts on `reopenedProject.options['proteomics.superseded_by']`. |
| A5 | `Project.name` (used for the `Proteomics-Review-<slug>-v<N>-<date>` convention) survives save + reopen. | PUB-09 + Pattern 1 | LOW — `Entity.name` is the canonical persistent identifier (`entities/entity.ts:40-41`). The live server's `recent projects` data shows project `name`/`friendlyName` are persistent across saves. |
| A6 | The volcano viewer config (axis labels, color binding, threshold lines) saved via `tv.getInfo()` survives the save → open cycle and renders identically on the reviewer side. | Pattern 1 + success criterion 3 + Phase 15-specific pitfall | MEDIUM — Phase 13 round-3 evidence (commit `e527d07ba1`) shows Filters viewer config is partially stripped. The volcano viewer may or may not have the same issue. Wave 0 test specifically asserts the volcano renders with the locked magenta/cyan/gray color binding and labeled axes on reopen. If it doesn't, plan adds a post-open `recomputeVolcano()` call as a defensive re-init step. |
| A7 | `User.friendlyName` (inherited from `Entity.friendlyName`) is the correct user-friendly source for `published_by`. | Claude's-discretion in CONTEXT.md | LOW — `Entity.friendlyName` is platform-standard; HitTriage's `hit-triage-app.ts:408` uses `grok.shell.user.friendlyName` for the same audit-log purpose. |
| A8 | The `dapi.projects.filter('name like "Proteomics-Review-${slug}-v%"')` smart-filter syntax accepts `like` predicates and matches multiple results. | D-04 republish detection + Open Question 6 | MEDIUM — js-api comments link to `https://datagrok.ai/help/datagrok/navigation/views/browse#entity-search` for the syntax; the exact predicate vocabulary (`like`, `contains`, `=`, `and`, `or`) is not enumerated in the code. Wave 0 test for `findPriorShare` is the verification surface — if `like` is wrong, fall back to `dapi.projects.list()` + client-side regex filter. |
| A9 | The slug regex `[A-Za-z0-9._-]+` produces names that `dapi.projects.save` accepts without rejection. | D-01 slug rules | LOW — the regex matches the documented entity name charset across the platform. Worst case: server rejects, the dialog surfaces an error and the publish does not partially-complete (the round-trip gate hasn't run yet). |
| A10 | If `assertPublishedShape` rolls back the new Project, the volcano viewer (`tv` reference held in memory) is still alive and the user can fix the error and retry. | Pitfall 1 mitigation / Pattern 1 | LOW — `dapi.projects.delete(project)` deletes a server-side entity; the in-memory `tv` and source `df` are unaffected. |
| A11 | The 7-column allowlist (PUB-02: Protein ID, Gene, log2FC, p-value, adj.p-value, sig flag, direction) is the COMPLETE set the published volcano needs. If Phase 13's Subcellular Location color binding or Phase 13's Comparison filter is in use on the source view, the trim must include those too. | Trim contract + Phase 15-specific pitfall | MEDIUM — REQUIREMENTS.md PUB-02 is exhaustive ("only Protein ID, Gene, log2FC, p-value, adj.p-value, sig flag, and direction"). Strict reading: no Subcellular Location, no Comparison. Audit-implication question for the planner: the published view will look different from the expert's view if those columns are bound on the expert's volcano. Discuss-phase D-domain pin says the audience is biologist-consumer who shouldn't need Subcellular Location depth, so strict trim is probably right — but planner should confirm. |
| A12 | The metadata column shape (multiple typed single-value columns prefixed `_meta_`) is preferred to JSON-in-one-column. | Pattern 3 + Claude's-discretion | LOW — column inspector renders typed columns legibly; the reopen-time recovery code reads them as typed values (no JSON.parse failure mode). The fallback (JSON-in-one-column) is captured in CONTEXT.md if column count gets unwieldy. 8 columns is not unwieldy. |
| A13 | Spectronaut Candidates source (which sets `proteomics.de_complete` at parse time) flows through `publishAnalysis` identically to non-Candidates source (which sets `de_complete` after DE pipeline). | D-05 + Wave 0 test fixtures | LOW — the only difference is which tags are set; `publishAnalysis` reads `proteomics.de_complete` (both shapes set it) and `proteomics.de_method` (Spectronaut Candidates sets `'spectronaut'`, others set `'limma'` etc. per CLAUDE.md tag table). The trim allowlist is the same. Wave 0 test explicitly covers both fixtures per Claude's-discretion in CONTEXT.md. |

**If a planner or operator wants to lower the risk on A1, A2, A3, A4, A6, A8, A11:** the Phase 15 plan's Wave 0 should include a one-shot exploration task (`publish-spike.ts`) that runs the full save→delete→find→open cycle against a tiny fixture and prints the actual shapes of `permissions.get`, post-rollback `tables.find`, `project.options` post-reopen, viewer round-trip, and smart-filter results. The spike's output feeds the locked decision for each assumption; the spike itself does not ship.

## Open Questions

These are decisions deferred to the per-task plan; they are not blockers.

1. **Smart-filter syntax for republish detection (`dapi.projects.filter`).**
   - What we know: `HttpDataSource.filter(w)` accepts a "smart filter" string per js-api/src/dapi.ts:325-332, linked to `https://datagrok.ai/help/datagrok/navigation/views/browse#entity-search`.
   - What's unclear: whether `name like "...-v%"` is the correct shape, or whether `name contains "Proteomics-Review-<slug>"` is preferred, and whether server-side tag filtering on `proteomics.published_target` is supported in addition to the name-pattern filter.
   - Recommendation: Wave 0 test runs both shapes against the live server with a fixture; the working shape becomes the locked syntax in `findPriorShare`. Fallback to `list()` + client-side filter is acceptable for small expected list sizes (<100 published projects per target).

2. **Slug uniqueness within a day for the same target.**
   - What we know: PUB-09 mandates `Proteomics-Review-<slug>-v<N>-<YYYY-MM-DD>` naming; D-04 mandates supersede-on-republish.
   - What's unclear: whether the slug+version combo is enough to disambiguate within a single day (republishing v2 twice on the same day produces `...-v2-2026-06-07` then `...-v3-2026-06-07` — fine), OR whether two separate analysts publishing the same target on the same day will collide (different version chains).
   - Recommendation: Use `findPriorShare(target, group).maxVersion + 1` as the version source. If two analysts publish concurrently, the second one's version increment is still correct as long as they both queried the prior. Add a save-failure retry with `+1` increment as a defensive layer.

3. **`childSpace` discovery when it already exists.**
   - What we know: `umbrellaClient.subspaceExists(name)` returns boolean; `umbrellaClient.addSubspace(name)` returns the created Space.
   - What's unclear: there's no `subspaceById(name)` method visible in js-api. Recovery path when subspace exists is: list children, filter by name. Slow O(N) for large umbrella.
   - Recommendation: Use `umbrellaClient.children.filter('Project', false).list()` and client-side filter by `friendlyName`. For Phase 15 the umbrella has ≤20 child Spaces typically (one per target); not a performance concern.

4. **Audit panel auto-dock on Project open.**
   - What we know: `@grok.decorators.panel` registers a context panel that fires on cell selection by `semType`.
   - What's unclear: whether there's a Project-open event or first-paint hook that auto-docks a panel without requiring a click. CONTEXT.md proposes `autostartImmediate`-style first-paint hook OR registration of a `meta.role: dashboard` function that runs on project open.
   - Recommendation: Two-step approach. Primary: `@grok.decorators.panel` with `semType = PROTEIN_ID` filter + `isPublished(df)` first-line guard — same shape as existing `uniprot-panel.ts`. The panel renders when the reviewer clicks any protein row, which is the natural first interaction. Optional enhancement: register an `@grok.decorators.func` with `meta.role: 'init'` or `autostartImmediate` that, on project open, calls `grok.shell.tv.dataFrame.currentRowIdx = 0` to trigger the panel automatically. Plan picks based on plan-time test of platform behavior.

5. **Mailto body templating.**
   - What we know: PUB-13 wording from CONTEXT.md: subject `"Re-run request: <published project name>"`, body `"Hi <sharer friendly name>, could you re-run with [different parameters]? Looking at <project name> (published <date>)."`.
   - What's unclear: how to populate the sharer's email address. `User.email` exists per `js-api/src/entities/user.ts:51`. But the sharer is the PUBLISHING user, whom the reviewer needs to email — the reviewer doesn't have direct access to the publishing user's email by default.
   - Recommendation: At publish time, store the sharer's email in `proteomics.published_by_email` tag (server-authoritative via `grok.shell.user.email`). The audit panel reads it. If `email` is null (some accounts), fall back to `mailto:?subject=...&body=...` (no recipient — opens user's mail client with empty To).

6. **`Project.options` vs `DataFrame.setTag` for the `superseded_by` pointer.**
   - What we know: `Project.options` is a `MapProxy` (entities/project.ts:26). DataFrames have `setTag/getTag` separate from this.
   - What's unclear: where to store the bidirectional supersede pointer. Two candidates: (a) on the inner `DataFrame.tags` (via `frozen.setTag('proteomics.superseded_by', newId)`) — discoverable when the DF is open in a view; (b) on the `Project.options` map (via `project.options['proteomics.superseded_by'] = newId`) — survives DF-content drift but requires fetching the Project entity.
   - Recommendation: Both. Set on the DF tag (so the volcano viewer code can read it without a Project fetch) AND on `project.options` (so an admin browsing the project list can see the chain without opening each one). Belt-and-braces — same philosophy as Pitfall 3.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| Datagrok server | Every step (`dapi.*`) | ✓ | `release/1.27.3` (live verified) | None — required |
| Datagrok js-api `Project.addChild`, `permissions.{grant,get}`, `spaces.{createRootSpace,id}`, `SpaceClient.addSubspace` | All publish steps | ✓ | `^1.25.0` declared; live `1.27.3` | None — required |
| Nested Spaces | D-03 primary shape | ✓ | Confirmed in live server response (namespace path `DDI:Target1:`) AND in js-api `addSubspace` API | D-03 fallback (Project nesting under flat Space) — NOT NEEDED per verification |
| Compute environment (R for limma/deqms) | Not used by Phase 15 | n/a | n/a | n/a — Phase 15 reads tags only, no R |
| `grok` CLI (for testing) | Wave 0 test development | ✓ | Available; live server reachable | Manual testing in browser if CLI unavailable |
| `@datagrok-libraries/test ^1.1.0` | `src/tests/publish-roundtrip.ts` | ✓ | Already installed | None — required |

**Missing dependencies with no fallback:** None.
**Missing dependencies with fallback:** None.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | `@datagrok-libraries/test ^1.1.0` (already installed) |
| Config file | None — tests registered via imports in `src/package-test.ts` |
| Quick run command | `grok test --category projects --test publish-roundtrip` |
| Full suite command | `grok test` (runs everything in `src/package-test.ts`) |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| PUB-01 | Publish menu item is registered + gated on `de_complete` | smoke | `grok test --test "PUB-01 menu visible"` | ❌ Wave 0 |
| PUB-02 | Trimmed snapshot has exactly 7 columns + Direction; no intensities | unit | `grok test --test "trim allowlist"` | ❌ Wave 0 |
| PUB-03 | Deep clone — source mutation does NOT leak (Pitfall 1) | unit | `grok test --test "frozen snapshot isolation"` | ❌ Wave 0 |
| PUB-04 | `proteomics.published_target` tag set; reachable post-reopen | integration | `grok test --test "target tag round-trip"` | ❌ Wave 0 |
| PUB-05 | View-only ACL on reviewer group; verify-and-rollback when Edit present (Pitfall 2) | integration | `grok test --test "verify-and-rollback"` | ❌ Wave 0 |
| PUB-06 | Volcano renders on open with stored thresholds | integration | `grok test --test "volcano round-trip"` | ❌ Wave 0 |
| PUB-07 | Audit panel shows DE method, thresholds, group names, target, date, sharer | integration | `grok test --test "audit panel content"` | ❌ Wave 0 |
| PUB-08 | Dialog captures target + group + note + confirmation | manual-only | n/a (UAT) | n/a |
| PUB-09 | Versioned/dated name survives reopen | integration | `grok test --test "project name round-trip"` | ❌ Wave 0 |
| PUB-10 | Republish creates new Project; supersede pointer set bidirectionally | integration | `grok test --test "supersede chain"` | ❌ Wave 0 |
| PUB-11 | Belt-and-braces metadata column survives round-trip | integration | `grok test --test "metadata column round-trip"` | ❌ Wave 0 |
| PUB-12 | Enrichment DF carried when present; cross-DF highlight re-wires | integration | `grok test --test "enrichment carry"` | ❌ Wave 0 |
| PUB-13 | Mailto link built with correct subject and body | unit | `grok test --test "mailto template"` | ❌ Wave 0 |

### Sampling Rate
- **Per task commit:** `grok test --category publishing` (only the publish round-trip suite — ~30 seconds with a fixture)
- **Per wave merge:** `grok test --category publishing,projects` (publish suite + Bio-pattern projects sanity)
- **Phase gate:** `grok test` (full suite green) before `/gsd:verify-work`

### Wave 0 Gaps

- [ ] `src/tests/publish-roundtrip.ts` — NEW, covers PUB-02, PUB-03, PUB-04, PUB-05, PUB-06, PUB-07, PUB-09, PUB-10, PUB-11, PUB-12, PUB-13. Two fixtures: synthetic demo (already used in `src/tests/analysis.ts`) + Spectronaut Candidates fixture (covers `de_complete` shortcut). Both fixtures round-tripped end-to-end with `assertPublishedShape` invoked.
- [ ] `src/publishing/assert-published-shape.ts` — NEW, single helper consumed by both `publishAnalysis` (defensive in-process check, Pitfall 3 gate) AND `publish-roundtrip.ts` test.
- [ ] (Optional spike) `src/tests/publish-spike.ts` — NEW, runs once during Wave 0 to verify Assumptions A1-A4 + A6 + A8 against the live server. Output: print actual shapes; do NOT ship.
- [ ] No framework install needed — `@datagrok-libraries/test` is already declared in package.json:22.

*Existing test infrastructure (`src/tests/volcano.ts`, `src/tests/parsers.ts`, etc.) is well-established; Phase 15 fits the same pattern with a new file plus an import line in `src/package-test.ts`.*

## Nested Spaces Verification

Per CONTEXT.md "D-03 PLAN-TIME RESEARCH FLAG (load-bearing)", confirm whether nested Datagrok Spaces are supported.

**Verdict: ✓ Nested Spaces are first-class.** Plan against the D-03 PRIMARY SHAPE.

**Evidence:**
1. **Live server data** [VERIFIED at 2026-06-07 against `release/1.27.3`]: `grok s raw GET /api/projects/recent` returns:
   ```
   Project  DDI:Target1:  7a66a5f0-...  ParentShortcutSpectronautHyeMix  target1-prot-2026-05-11
   ```
   The `namespace` column shows `DDI:Target1:` — a Project nested inside Space `Target1`, which itself is nested inside Space `DDI`. This is exactly the D-03 primary shape (umbrella `DDI` → child `Target1` → Project for one target's analysis).
2. **`js-api/src/dapi.ts:805-866`** (verified by direct read): `SpacesDataSource` + `SpaceClient` + `SpaceChildrenClient` + `SpaceFilesClient` — full API surface for create, addSubspace, addEntity, exists, files. Comments explicitly describe hierarchy: *"Adds a child space (subspace) under this space. The subspace will have its own storage with a namespaced path derived from the parent hierarchy (e.g., `'ParentSpace:ChildSpace:'`)"*.
3. **`packages/ApiSamples/scripts/dapi/spaces.js:38-48`** (verified): live demo of `spaceClient.addSubspace('ChildSpace')` and `spaceClient.files.createDirectory(subspaceName)` (directory creation auto-creates subspace).
4. **Cleanup behavior** (per `spaces.js:88`): `await grok.dapi.spaces.delete(space)` — *"cascades to subspaces and files"*. So rollback at the Space level cleanly removes child Spaces and contained Projects.

**Plan implication:** Drop the D-03 fallback section from the plan. Use the primary shape directly. Recommend the planner explicitly note "D-03 fallback NOT NEEDED — nested Spaces verified" in the PLAN.md so a future reader doesn't re-litigate.

## Pitfall 3 Round-Trip Enumeration

Per CONTEXT.md "Pitfall 3 round-trip enumeration (load-bearing)", confirm which `proteomics.*` tags + `Proteomics-*` semTypes + `df.name` survive `DG.Project.save` followed by `.open()`.

**Status:** Enumeration script NOT executed at research time (requires live save+reopen against the active worktree, which would create test entities in the live server). Recommendation: planner generates `src/tests/publish-spike.ts` in Wave 0 — this is a one-shot exploration test that runs once, prints the shapes, and informs the locked decisions. Below is the known evidence and the expected enumeration shape.

**Known-stripped categories (HIGH confidence per Phase 13 evidence):**
- Layout viewer config — Phase 13 commit `e527d07ba1` proves the platform serializer strips `look.filters[]` and `look.columnNames` from a docked Filters viewer regardless of input shape. The trim contract MUST NOT depend on Filters viewer config surviving; the published view re-docks the Filters viewer post-open (planner gates on PUB-06 + the carry-forward from Phase 14 D-05).

**Likely-surviving categories (MEDIUM confidence):**
- Project name + friendlyName — `Entity` base provides these as durable per `entities/entity.ts:36-41` and live server data shows them stable across renames. [VERIFIED: live data]
- DataFrame contents (columns + cell values) — primary purpose of `tables.save`; trivially survives. [VERIFIED: Bio test pattern]
- Column data types — required for table reconstruction; survives. [ASSUMED: not explicitly verified, but trivially true since the test framework reads columns by name post-reopen]
- `df.name` (set on the trimmed clone before save) — `DataFrame.name` is persisted as part of TableInfo; survives. [ASSUMED based on Bio test which reopens by id and reads `tv.dataFrame.name`]

**Unknown-survival categories (the spike will enumerate):**
- `proteomics.*` DataFrame tags (`proteomics.de_method`, `proteomics.groups`, `proteomics.published`, etc.) — Phase 13 evidence is about viewer-config tags; pipeline tags may survive or may not. The belt-and-braces metadata column eliminates the dependency, but the spike confirms whether the tags ALSO survive (in which case the column is defensive but not critical).
- `Proteomics-*` semTypes on individual columns — these drive the UniProt panel auto-fire (SEMTYPE.PROTEIN_ID), the volcano binding, the audit panel filter. The `detectors.js` mirror MIGHT auto-reapply on reopen (CLAUDE.md mentions detectors load before bundle); the spike confirms whether detectors re-fire post-reopen.
- `Project.options` MapProxy entries — used in plan for `superseded_by` (Open Question 6). The spike confirms readability post-reopen.
- Viewer config: volcano color binding, axis labels, threshold lines, top-N labels (Phase 14 D-03 `setTopNLabels`) — load-bearing per Phase 15-specific pitfall.

**Spike implementation outline (for Wave 0):**
```typescript
// src/tests/publish-spike.ts
import {category, test} from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

category('publish-spike', () => {
  test('save+open round-trip — enumerate survivors', async () => {
    const df = makeMinimalDeFixture();   // 10 proteins, log2FC, p-value, adj.p, sig
    df.setTag('proteomics.de_method', 'limma');
    df.setTag('proteomics.published_target', 'MYH7');
    df.col('Protein ID')!.semType = 'Proteomics-ProteinId';
    df.name = 'spike-source';

    const view = grok.shell.addTableView(df);
    const volcano = view.scatterPlot({x: 'log2FC', y: 'adj.p-value'});

    const proj = DG.Project.create();
    proj.name = 'spike-test';
    proj.options['custom_option'] = 'hello';
    proj.addChild(df.getTableInfo());
    proj.addChild(view.getInfo());
    await grok.dapi.tables.uploadDataFrame(df);
    await grok.dapi.tables.save(df.getTableInfo());
    await grok.dapi.views.save(view.getInfo());
    await grok.dapi.projects.save(proj);

    const projId = proj.id;
    grok.shell.closeAll();

    const reopened = await grok.dapi.projects.find(projId);
    await reopened.open();
    const reDf = grok.shell.tv.dataFrame;

    // Enumerate and PRINT (do not assert)
    console.log('df.name survived:', reDf.name);
    console.log('de_method tag:', reDf.getTag('proteomics.de_method'));
    console.log('target tag:', reDf.getTag('proteomics.published_target'));
    console.log('PROTEIN_ID semType:', reDf.col('Protein ID')?.semType);
    console.log('options.custom_option:', reopened.options['custom_option']);
    console.log('viewers:', grok.shell.tv.viewers.map((v) => v.type));

    // Cleanup
    await grok.dapi.projects.delete(reopened);
  });
});
```
This spike runs ONCE, prints the shapes, the planner reads the output and locks the trim+belt-and-braces contract accordingly. The spike does not ship in the final phase.

## Sources

### Primary (HIGH confidence)
- `js-api/src/dapi.ts` at HEAD `21d37776ee` — full Dapi surface: `projects.save/find/delete/filter`, `permissions.grant/get/check/revoke`, `groups.list/filter`, `spaces.createRootSpace/rootSpaceExists/id`, `SpaceClient.addSubspace/subspaceExists/addEntity/children/files`, `users.current`, `HttpDataSource.filter/list/find/delete/save`
- `js-api/src/entities/project.ts` at HEAD — `Project.create`, `addChild`, `addLink`, `description`, `options` (MapProxy), `open(options)`, `path`, `isSpace`
- `js-api/src/entities/entity.ts` at HEAD — `Entity.friendlyName`, `name`, `id`, `nqName`
- `js-api/src/entities/user.ts` at HEAD — `User.current`, `friendlyName`, `email`, `group`
- `packages/ApiSamples/scripts/dapi/projects.js` — canonical 3-line publish shape
- `packages/ApiSamples/scripts/dapi/spaces.js` — nested-Space + addSubspace + addEntity + delete-cascades
- `packages/ApiSamples/scripts/dapi/layouts-and-permissions.js` — `permissions.grant(entity, group, edit:bool)` + `permissions.get(entity)` pattern
- `packages/Bio/src/tests/projects-tests.ts:26-47` — multi-`addChild` (tableInfo + layoutInfo), save→find→open round-trip, viewer presence assertion pattern
- `packages/HitTriage/src/app/hit-triage-app.ts:370-429` — AppData persistence + `grok.shell.user.friendlyName` audit usage + `_friendlyName` precedence for repeat-save author preservation
- `packages/Proteomics/src/viewers/volcano.ts` — Phase 14 D-04 color lock implementation (`DIRECTION_COLORS_BASE`), `VOLCANO_METRIC_TAG`, `ensureNegLog10Column`, `ensureDirectionColumn` — load-bearing for Phase 15 round-trip carry
- `packages/Proteomics/src/viewers/enrichment-viewers.ts:9` — module-level subscription pattern for cross-DF protein-highlight (D-05 reuse)
- `packages/Proteomics/src/viewers/qc-computations.ts:27-32` — `ensureFreshFloat` idempotency pattern (carries to trim helper)
- `packages/Proteomics/src/utils/proteomics-types.ts` — verified SEMTYPE set; no new types added in Phase 15
- `packages/Proteomics/src/utils/column-detection.ts` — `findColumn` + `findProteomicsColumns` (trim allowlist resolves through these)
- `packages/Proteomics/CLAUDE.md` — pipeline tag conventions, `findColumn`/`SEMTYPE` requirements, clone-for-isolation precedent, function-naming prefixes
- `packages/Proteomics/.planning/codebase/ARCHITECTURE.md` §"In-place Mutation" + §"Anti-Patterns" §"Mutation of the In-Place DataFrame Without Idempotency Guard"
- `packages/Proteomics/.planning/codebase/STRUCTURE.md` §"Where to Add New Code"
- `packages/Proteomics/.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md` — locked decisions D-01..D-05 + Claude's discretion
- `packages/Proteomics/.planning/REQUIREMENTS.md` PUB-01..PUB-13
- `packages/Proteomics/.planning/research/{SUMMARY,STACK,ARCHITECTURE,PITFALLS}.md`
- `packages/Proteomics/.planning/milestones/v1.3-phases/14-ck-omics-analyst-experience-enhancements/14-CONTEXT.md` — Phase 14 D-04 / D-05 carry-forward
- `packages/Proteomics/.planning/milestones/v1.3-phases/13-ck-omics-volcano-and-enrichment-parity/13-CONTEXT.md` — Phase 13 evidence anchor for Pitfall 3
- `tools/GROK_S.md` lines 23, 144-152, 253-265 — `grok s shares add/list`, raw API, namespace path conventions
- Live Datagrok server (`release/1.27.3`, `grok s raw GET /api/projects/recent` + `/api/spaces`) — verified at 2026-06-07 against the configured local instance; confirms nested Spaces and project namespace path convention

### Secondary (MEDIUM confidence)
- `packages/HitTriage/src/app/dialogs/permissions-dialog.ts` — referenced from `hit-triage-app.ts:18`; not read in detail but identified as the precedent for a permission-editing UI (Phase 15 does NOT need this — D-02 picks group from `dapi.groups.list()` and a single `permissions.grant` call)
- Saved memory `feedback_keep_workaround_capture_future` — belt-and-braces metadata column is a workaround; track cleanup as a future-action todo when platform serializer is fixed

### Tertiary (LOW confidence — flagged in Assumptions Log)
- Smart-filter expression syntax for `dapi.projects.filter` (Open Question 1 / Assumption A8) — js-api comments link to docs at `https://datagrok.ai/help/datagrok/navigation/views/browse#entity-search` but the predicate vocabulary is not in the code
- Cascade behavior of `dapi.projects.delete` (Assumption A3)
- Permission inheritance behavior when moving a Project INTO a Space that already has grants (Assumption A2)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — every API verified at js-api HEAD AND against live server `1.27.3`
- Architecture: HIGH — Phase 15 reuses well-established v1.3 patterns; new code goes into `src/publishing/` per ARCHITECTURE.md contract
- Pitfalls: HIGH for the four well-documented ones (1, 2, 3, 4, 14) + one new Phase 15-specific pitfall (volcano carry-forward) — all anchored in repo evidence or live server behavior
- Validation: HIGH — `@datagrok-libraries/test` already installed, existing test file shape is well-precedented in `src/tests/volcano.ts`; round-trip pattern is directly transposable from Bio
- D-03 nested Spaces: HIGH — verified in BOTH js-api code AND live server data
- Pitfall 3 tag-survival enumeration: MEDIUM — known categories documented; final enumeration deferred to Wave 0 spike script

**Research date:** 2026-06-07
**Valid until:** 2026-07-07 (30 days for the platform surface — js-api and live server `1.27.3` are stable; Phase 15 plan should land within this window or re-verify if the platform releases a minor)

---
*Phase: 15-read-only-publishing-foundation*
*Researched: 2026-06-07*

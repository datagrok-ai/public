# Compute2 Artifact Alignment PoC — Implementation Plan

Implements the design in `artifact-alignment.html`, narrowed to a single artifact type:
**saved Compute2 workflow runs**. The doc stays general; every compute2-specific decision
lives here.

## Scope and decisions (from the design discussion)

- **Versioning**: latest-approved-stays-live (as now specified in the doc) — a pending or
  rejected update never affects what the audience sees.
- **Republish key**: `(program, study, name)`.
- **Curation**: mutable post-publish property schema — `tags` (string array), `path` (string).
- **FuncCall permissions are ignored in the PoC**: the publish flow contains the grant/revoke
  calls **commented out**, marking exactly where core enforcement slots in once FuncCall becomes
  a full entity. Everything else about the frozen copy (TableInfo grants, group scoping,
  audit) works as designed.
- **No client-side hiding hacks** (no soft-deleted TableInfos, no blob-storage detour) — they
  regress permission handling. Hiding nested entities is a **core ask** (see the end).
- **Gate OFF** for workflow runs → publish auto-approves; the reject/approve path is still
  implemented and exercisable manually (needed for tests and for gated types later).
- Out of scope: master-data sync, bulk publish, save/share dialog injection (context-menu /
  ribbon action instead — the share dialog is closed to packages), true service identity
  (emulated, see Phase 2), all other artifact types.

## What a compute2 artifact is (the object being published)

One saved workflow run = a graph:

- **meta FuncCall** of the pipeline's nqName — `options[PIPELINE_CONFIG]` holds the serialized
  tree with per-step `funcCallId` references
  (`reactive-tree-driver/src/runtime/funccall-utils.ts:120`, `PipelineInstance.ts:106`);
- **one FuncCall per step** (`StateTree.save()` → `saveFuncCall` → `historyUtils.saveRun`);
- **one TableInfo per DataFrame input/output of every step** plus consistency-restriction
  DataFrames (`_DG_CONSISTENCY_DF_REF_` refs inside `options[INPUT_RESTRICTIONS]`);
- file-input blobs written via `grok.dapi.files.write`.

`artifact_id` in the alignment row = the **frozen meta call id**; `source_artifact_id` = the
author's meta call id.

## Deliverable shape

New package **`ArtifactAlignment`** (working title), depending on
`@datagrok-libraries/compute-utils`:

```
ArtifactAlignment/
  src/
    schema/        domain schema deploy + group provisioning
    service/       publish, approve/reject, republish, drift check, graph clone
    app/           catalog Browse app (program landing, facets, history pane)
    compute2/      publish dialog + open-readonly glue
  package.ts       #app + registered functions
```

Compute2 itself gets only a thin, **optional** hook (see Phase 3): the publish action appears
only when the `ArtifactAlignment` package is installed **and** a Compute2 package-level
integration setting is enabled — no hard dependency in either direction.

## Phase 1 — Domain schema and groups

1. Registries: `program` (with viewers/contributors/approvers group refs), `compound`,
   `program_compound`, `study` — manual admin path only (no sync).
2. `alignment` table per the doc's model: identity columns immutable;
   `artifact_type` fixed to `workflow-run`; unique `(publication_id, status)` among live rows;
   `study` part of the republish key (null study is a distinct key value).
3. `approval` property schema — `status` default `'pending'`, Edit → approvers group,
   View broad.
4. `curation` property schema — `tags: string[]`, `path: string`, Edit → contributors group,
   View broad.
5. `alignment_history` — businessKey `(publication_id, revision)`, full column copy +
   `superseded_on`, `original_row_id`.
6. Group provisioning function per program; write-coupling accepted per the doc (contributors
   get Edit on the program row).

All mechanics here were verified on a live deployment per the doc (choices, delegate security,
property-schema grants via the domains REST route, M:N, businessKey idempotence). The only new
constraint shape is the composite live-row unique key — verify it first (same mechanism class
as the verified single-column unique; a one-hour fixture probe).

## Phase 2 — Publish service (deep clone of the run graph)

`publishWorkflowRun(sourceMetaCallId, {program, study, name, workstream, compounds, description})`:

1. `loadInstanceState(sourceMetaCallId)` → serialized config + source meta call.
2. Mint the frozen meta call id up front. For each `funcCallId` in the config:
   `historyUtils.loadRun` (materializes DFs and files), `fc.newId()`, stamp
   `options.parentCallId = <frozen meta call id>`, save through a publish variant of `saveRun`
   that grants each uploaded TableInfo **View to the program viewers group** (not All users).
3. Re-upload consistency DataFrames the same way; rewrite the `_DG_CONSISTENCY_DF_REF_` ids.
4. Rewrite `funcCallId` refs in the config; `newId()` the meta call; stamp
   `options.publicationId / revision / frozen='true'`; save.
5. **Commented permission calls** on every frozen funccall:
   ```ts
   // Requires core: FuncCall ACL + immutability by ownership.
   // await grok.dapi.permissions.grant(frozenCall, viewersGroup, false);
   // await catalog.transferOwnership(frozenCall, serviceUser);
   ```
6. Alignment row: resolve publication by `(program, study, name)`. First publish → mint
   `publication_id`, revision 1. Republish → archive any stale live in-review row (idempotent
   history insert), insert revision+1 **alongside** the live approved row, copy curation
   forward from the latest version.
7. Approval step (auto, gate off): archive the previously approved row → soft-delete it → set
   the new row `approved`. Every partial-failure window is idempotently repairable; the drift
   check (below) scans for interrupted flips and for publications with no live row.

**Service identity is emulated**: the PoC runs as the publishing user, so frozen entities are
author-owned. Immutability rests on RTD's readonly load path plus the drift check until core
provides funccall ACL/ownership transfer. This is the PoC's biggest declared gap.

Also in this phase: `driftCheck()` — grants mismatch (intended vs. actual on frozen
TableInfos), live-row invariants, frozen-graph integrity (every `funcCallId`/table ref
resolves).

## Phase 3 — Compute2 integration (optional, double-gated)

The integration is off unless **both** hold:

- the `ArtifactAlignment` package is installed —
  `ArtifactAlignment:publishWorkflowRun` resolves via `DG.Func.find`;
- Compute2's package-level setting `enableArtifactPublishing` (a `properties` entry in
  Compute2's `package.json`, default `false`, editable in the standard package settings UI) is
  turned on.

Both checks run once at action-registration time; either one failing means Compute2 renders
exactly as today (no ribbon item, no dialog code loaded — the glue lives behind a dynamic
import).

1. **Publish action** in TreeWizard (ribbon / global action, shown when the workflow has a
   saved state): dialog with Name / Program / Study / Molecules / Workstream / Description,
   pre-filled from the run's metadata; republish hint per the doc ("publishes version N; the
   current version stays visible until this one is approved"). Compounds are a plain
   multi-select in the PoC (no semantic-type auto-suggestion).
2. **Open path**: the catalog opens the frozen meta call through the existing readonly load
   (`loadFuncCall(id, isReadonly=true)` / Compute2 app URL). Verify readonly is forced for the
   whole tree when loading a frozen instance — add a driver-level `forceReadonly` load flag if
   the serialized state doesn't already carry it.
3. Optional (defer if tight): "Clone to my workspace" from a frozen artifact — the same clone
   routine without the alignment insert, producing a mutable working copy.

## Phase 4 — Catalog UI (Browse app)

Minimal slice of the doc's main page:

- Programs tree → program landing: approved-only default, grouped by workstream; "in
  re-review" facet off by default (author + approvers see their rows).
- Facets from the shipped domain filter UI: status, study, workstream, molecule (M:N), type,
  **tags**, **path**.
- Item click → Compute2 readonly open. Context panel: alignment metadata, curation fields
  editable in place (contributors), History pane from `alignment_history` with per-version
  open.
- Routes by `publication_id`.
- Deferred: "Open from program" picker dialog, endorsement ranking, global-search providers.

## Phase 5 — Tests (headless-first)

Package tests at the funcCall/API level (no UI driving), with server-state cleanup in
`tearDownAll`:

- publish → approved-slice query returns exactly one row; frozen graph loads readonly and
  structurally equals the source (serialized-state snapshot, ids masked);
- republish same key → revision+1 pending while the previous row stays approved-live;
  approve → old row archived (history businessKey present), audience query flips atomically
  enough (repairable window covered by a drift-check test);
- reject → audience query unchanged; resubmit → rejected row archived;
- widened key: same name in a different study → distinct `publication_id`;
- curation: contributor edits tags/path; values copied forward on republish; facet query by
  `tags`/`path` returns the row; non-contributor write → typed denial;
- gating: with the package absent or `enableArtifactPublishing` off, Compute2 registers no
  publish action; with both on, the action appears;
- grants: a viewers-group member opens the frozen run and all its tables; a non-member is
  denied on the tables (funccall shells load regardless — asserted and documented as the known
  core gap);
- idempotence: duplicate history insert no-op; duplicate live insert typed unique violation;
  drift check repairs an artificially interrupted approve flip.

## Known gaps accepted in the PoC → core asks

| Gap | PoC behavior | Core ask |
|---|---|---|
| FuncCall is not a full entity (no author column, tags, ACL, immutability) | readonly-by-UI, public-by-id; commented grant calls mark the slots | Promote FuncCall: author/createdOn, tags, per-entity ACL, immutability by ownership |
| Nested entities pollute the author's My stuff (run TableInfos) | unchanged for working saves; published copies are group-granted and (once service identity exists) service-owned | Nested/hidden entity support: entities ownable by a parent (funccall/publication) and excluded from galleries, recents, My stuff |
| `uploadDataFrame` always creates a TableInfo | accepted | Expose the Dart `saveInfo: false` parameter in js-api |
| No service identity from a client plugin | publish runs as the author; ownership transfer TODO-marked | Service-identity execution or ownership-transfer API |
| Share dialog closed to packages | ribbon/context action | Entity dialog extension points (already in the doc's blockers) |
| Working saves still grant All users on tables | unchanged | Follows from nested-entity hiding + funccall ACL |

## Order of work

Phase 1 (with the composite-key probe first) → Phase 2 + its tests → Phase 3 (optional,
double-gated) → Phase 4 → Phase 5 completes alongside each phase, not at the end. Phases 1–2
and 4 are fully exercisable without the Compute2 integration (publish via a registered
function / the catalog app), so Phase 3 can land last or be deferred without blocking the PoC.

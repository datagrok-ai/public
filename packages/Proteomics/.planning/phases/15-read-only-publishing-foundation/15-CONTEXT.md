# Phase 15: Read-Only Publishing Foundation - Context

**Gathered:** 2026-06-07
**Status:** Ready for planning

<domain>
## Phase Boundary

A trimmed, frozen, target-keyed, view-only-shared `DG.Project` snapshot of a DE-complete proteomics analysis that **consuming biologists** (NOT fellow proteomics practitioners) can open via a direct link, see a first-paint volcano + audit context panel, and survive the Datagrok platform serializer's known tag-stripping behavior. Republish creates a new versioned Project; the previous version gets a `proteomics.superseded_by` pointer and stays available for audit. View-only ACL is enforced post-grant via `dapi.permissions.get` verification with hard rollback if Edit slips in via Space inheritance.

Scope is fixed by REQUIREMENTS.md PUB-01..PUB-13 (P1 = PUB-01..PUB-11; P2 = PUB-12 enrichment carry, PUB-13 mailto). This discussion clarifies HOW to implement what's locked, not WHETHER to add new capabilities. v1.4-cross-cutting capabilities (SPC tracking, campaign data model, cross-run comparison) belong in Phases 16/17/18 and MUST NOT be pulled into Phase 15.

**Audience pin (load-bearing):** "Reviewer" means consuming biologist throughout this phase. Drives the jargon audit (Pitfall 14), the trim allowlist's rationale (no raw intensities — a biologist won't spot-check at intensity level), the audit-context-panel framing (narrate what an expert would assume known), the "Request re-run with different parameters" mailto wording (PUB-13), and the publish-destination decision (D-03 below — biologists land via link, not by Space-browsing).

</domain>

<decisions>
## Implementation Decisions

### Target identifier shape (PUB-04, PUB-09)
- **D-01:** **Freeform string + slug sanitization for the project name; raw label preserved in tag + metadata column.** User types any string in the share dialog (e.g. `MYH7-DMD`, `Cytokinetics atrophy panel`). Publish helper slugifies to `[A-Za-z0-9._-]+` for the project-name slot in `Proteomics-Review-<slug>-v<N>-<YYYY-MM-DD>` (PUB-09). The raw user string is stored in BOTH `proteomics.published_target` tag AND a single-row metadata column (belt-and-braces, Pitfall 3 mitigation). No registry, no taxonomy, no constrained list — first-class target taxonomy is explicitly deferred to v1.5+ per REQUIREMENTS.md.
  - Slug rules (Claude's discretion to tune): preserve case; replace whitespace + special chars with `-`; collapse repeats; drop trailing punctuation; cap at ~64 chars to keep project names readable.

### Reviewer group sourcing (PUB-05, PUB-08)
- **D-02:** **`DG.ChoiceInput` populated from `grok.dapi.groups.list()`, filtered to groups the publishing user can administer.** User picks from a real list of existing groups — no typo path to a non-existent group, no dead-end share attempts. Matches platform muscle memory. Membership management stays upstream of this package (REQUIREMENTS.md "Out of Scope" — permission grants tied to a group, group lifecycle is admin-managed).

### Publish destination + inheritance enforcement (PUB-05; Pitfall 2 mitigation)
- **D-03:** **Nested Datagrok Spaces.** Umbrella Space `Proteomics-Reviews` contains per-target child Spaces `Proteomics-Review-<slug>`; each child Space holds the versioned Projects for that target (v1, v2, ...). Reviewer-group View grant lands at the **per-target child Space** level so all versions for that target inherit the grant uniformly. Post-grant verification gate is non-negotiable: `await grok.dapi.permissions.get(publishedProject)` MUST show the reviewer group's privilege list `=== ['view']`. If `Edit` / `Share` / `Delete` appears (Pitfall 2 — Space-inheritance leak), the publish helper THROWS, the new Project is rolled back (deleted), and the dialog surfaces `"Reviewer group already has elevated access via Space inheritance — publish aborted; ask an admin to scope the umbrella Space's permissions"`. Reviewer audience model (D-domain pin) is biologist-via-link, so the per-target Space is for admin audit / future-target-search, NOT a biologist discovery surface — biologists receive the project link directly from the expert.
  - **PLAN-TIME RESEARCH FLAG (load-bearing):** Confirm Datagrok supports **nested Spaces** as a first-class primitive (parent Space contains child Spaces, with permission inheritance flowing parent→child). If not supported, the fallback shape is Project-nesting under one flat Space: umbrella Space `Proteomics-Reviews` (flat) + per-target parent Project `Proteomics-Review-<slug>` (acts as a container) + child version Projects nested inside. Researcher MUST verify this against the live platform before the planner commits to the nested-Space contract; both candidate shapes are diagrammed below.

```
D-03 PRIMARY SHAPE (nested Spaces):
Space: Proteomics-Reviews              (umbrella)
  Space: Proteomics-Review-MYH7        (child, per target — View grant lands HERE)
    Project: Proteomics-Review-MYH7-v1-2026-06-07
    Project: Proteomics-Review-MYH7-v2-2026-06-14
  Space: Proteomics-Review-DMD
    Project: Proteomics-Review-DMD-v1-2026-06-08

D-03 FALLBACK SHAPE (if nested Spaces unsupported — Project nesting):
Space: Proteomics-Reviews              (flat, single umbrella)
  Project: Proteomics-Review-MYH7      (parent / container)
    Project: Proteomics-Review-MYH7-v1-2026-06-07
    Project: Proteomics-Review-MYH7-v2-2026-06-14
  Project: Proteomics-Review-DMD
    Project: Proteomics-Review-DMD-v1-2026-06-08
```

### Republish UX (PUB-10, PUB-08)
- **D-04:** **One menu entry `Proteomics → Share → Share Analysis for Review...` with detect-and-prefill behavior.** On dialog open, helper queries `grok.dapi.projects.filter` for matching `proteomics.published_target` + reviewer group (the active DF's last-share, if any, takes precedence). If a prior share exists:
  - Pre-fill target / reviewer group / note from the most recent matching share.
  - Render a banner at the top of the dialog: `"⚠ This will publish as v<N+1> and supersede '<prior project name>'"`.
  - On OK: helper sets `proteomics.superseded_by = <new project id>` on the prior version (and `proteomics.supersedes = <prior id>` on the new one — bidirectional audit trail). Prior version is NOT deleted (Pitfall 4 — soft pointer, never destructive).
  - Reviewer's bookmark to v1 still resolves; the reviewer-side audit panel on v1 surfaces a "Newer version available: [link]" note (Pitfall 14 UX hygiene — `published-analysis-panel.ts` checks for incoming `superseded_by` and renders the link).
  - No separate "Update Share" menu item — keeps menu surface flat (no Proteomics submenu today pairs create/update entries).

### Enrichment carry (PUB-12, P2)
- **D-05:** **Two DataFrames in the published Project.** Trimmed protein DF (the 7-column allowlist from PUB-02) + trimmed enrichment DF (term / source / p / adj.p / intersection-genes columns — Claude picks the exact allowlist from existing v1.2 enrichment columns). Cross-DF protein-highlight wiring re-establishes on Project reopen via the existing v1.2 `enrichDf.onCurrentRowChanged` subscription pattern (`src/viewers/enrichment-visualization.ts` precedent). This is **explicitly NOT a new `CampaignSelectionBus`** — that abstraction is Phase 18's territory (Pitfall 13 — campaigns are N-producer/N-consumer graphs; Phase 15 is single-producer/single-consumer and the v1.2 pattern fits).
  - When source DF has no enrichment results (no `proteomics.enrichment` tag or no enrichment DF in the shell), publish ships the protein DF only — the enrichment carry is opportunistic. PUB-12 is P2 (differentiator), so a no-enrichment source must still publish cleanly.
  - Round-trip test (Pitfall 3) covers BOTH the no-enrichment and with-enrichment shapes.

### Claude's Discretion
- **Direction column representation (Up / Down / NS).** String or numeric (±1/0/-1); both can drive the volcano color binding. Pick whichever the v1.3 volcano already consumes; if neither exists pre-publish, computed-at-publish-time string is simpler.
- **Sharer's friendly name source (PUB-07).** `grok.user.current().friendlyName` or equivalent; check what `dapi.users` returns and pick the user-friendly field. Stored in `proteomics.published_by` tag AND metadata column (belt-and-braces).
- **Single-row metadata column shape (PUB-11, Pitfall 3).** Multiple typed one-row columns (one per critical tag: `published_target`, `published_at`, `published_by`, `published_de_method`, `published_fc_threshold`, `published_p_threshold`, `published_version`, `published_id`) — column inspector renders them legibly and the reopen-time recovery code reads them as typed values. JSON-string-in-one-column is the fallback if the column count gets unwieldy.
- **Audit context panel triggering (PUB-07).** Register `published-analysis-panel.ts` as `@grok.decorators.panel` with `semType = PROTEIN_ID` filter + an `isPublished(df)` first-line check that early-returns `null` on non-published DFs. ALSO add a one-time `autostartImmediate` or first-paint hook that auto-docks the panel on Project open when `df.getTag('proteomics.published') === 'true'` — biologist shouldn't have to click around to find the audit context.
- **Round-trip test fixture choice (Pitfall 3).** Synthetic demo dataset already used in v1.0 tests (small, deterministic) + one Spectronaut Candidates fixture (covers the `de_complete` shortcut path where DE is pre-computed). Both shapes round-tripped end-to-end in `src/tests/publish-roundtrip.ts`.
- **Slug sanitization rules (D-01).** See D-01 notes; planner picks exact charset + cap.
- **Mailto body shape (PUB-13, P2).** Subject: `"Re-run request: <published project name>"`. Body: `"Hi <sharer friendly name>, could you re-run with [different parameters]? Looking at <project name> (published <date>)."` — pre-filled, biologist edits before sending. Plain `mailto:` URL — no email service integration.
- **`isPublished(df)` helper signature.** Single boolean tag check on `proteomics.published === 'true'` (set unconditionally by the publish helper alongside `published_at`). Used by the panel, by re-publish detection, and by any future code that wants to gate behavior on "this is a published clone, not the live source".

### Folded Todos
- **`2026-05-11-share-analysis-read-only-with-biologics-team-filed-by-target.md`** — explicit source todo for the whole milestone story; "filed by target" = PUB-04 target-keyed filing; "read-only" = PUB-05 view-only ACL + Pitfall 2 verify-and-rollback; "biologics team" = the biologist-consumer audience pin in the domain section. Captured here so the planner sees the original use-case wording, not just the abstracted REQ IDs.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase 15 requirements + scope (LOCKED)
- `.planning/REQUIREMENTS.md` PUB-01..PUB-13 — the locked requirements for this phase. P1 = PUB-01..PUB-11; P2 = PUB-12 (enrichment carry, D-05 above) + PUB-13 (mailto, Claude's discretion above). Do NOT discuss or replan requirement WHATs — they are fixed.
- `.planning/ROADMAP.md` Phase 15 entry — phase goal + 5 numbered success criteria. The round-trip test (`assertPublishedShape`) is the load-bearing gate per success criterion 3.
- `.planning/PROJECT.md` v1.4 Cross-Team Review section — milestone-level positioning + "biologist consumers" audience phrasing (drives D-domain pin).

### Research synthesis (load-bearing for plan-time research)
- `.planning/research/SUMMARY.md` §"Phase 15: Read-Only Publishing Foundation" + §"Research Flags" — four-phase ordering rationale and the explicit Phase 15 plan-time research task ("which `proteomics.*` tags + `Proteomics-*` semTypes + `df.name` survive `DG.Project` save/reopen"). Researcher MUST execute the throw-away round-trip enumeration script as the first plan-time deliverable.
- `.planning/research/PITFALLS.md` §"Pitfall 1" (stale-snapshot leak) / §"Pitfall 2" (Space-inheritance Edit slip — D-03 mitigation lives here) / §"Pitfall 3" (tag-stripping on Project serialization — D-04/D-05 belt-and-braces lives here) / §"Pitfall 4" (versioning ambiguity — D-04 supersede chain lives here) / §"Pitfall 14" (Cytokinetics demo audience contains biologists — domain pin + Claude's-discretion UX hygiene lives here). All five pitfall mitigations are baked into the decisions above; planner verifies each in the plan's task structure.
- `.planning/research/ARCHITECTURE.md` §"src/publishing/ NEW directory" — file layout contract for the new sibling dir (`publish-state.ts`, `trim-dataframe.ts`, `publish-project.ts`, `share-dialog.ts`); 5 new tags namespace (`proteomics.published`, `published_at`, `published_by`, `published_target`, `published_audit`) confirmed non-colliding via grep.
- `.planning/research/STACK.md` §"`DG.Project` + `grok.dapi.permissions.grant`" — canonical 3-line shape verified against `packages/ApiSamples/scripts/dapi/projects.js` and `packages/Bio/src/tests/projects-tests.ts:26`.

### Phase 13/14 carry-forward (do NOT regress)
- `.planning/milestones/v1.3-phases/14-ck-omics-analyst-experience-enhancements/14-CONTEXT.md` — Phase 14 D-04 (magenta/cyan/gray volcano color lock) + D-05 (Filters viewer scoping pattern) MUST survive into the published volcano. The published view should look identical to the expert's view (modulo trim) so the biologist sees the same colors the expert showed them in person.
- `.planning/milestones/v1.3-phases/13-ck-omics-volcano-and-enrichment-parity/13-CONTEXT.md` — Phase 13 D-02 (userDataStorage cross-session cache pattern) is the precedent for any per-publishing-user state Phase 15 needs. Commit `e527d07ba1` is the serializer-strip evidence anchor for Pitfall 3.
- `.planning/codebase/ARCHITECTURE.md` §"In-place Mutation" + §"Cross-DataFrame Selection" + §"Anti-Patterns" — confirms why deep-clone-first (Pitfall 1) is non-negotiable for publish; confirms the v1.2 enrichment subscription pattern D-05 reuses.
- `.planning/codebase/STRUCTURE.md` §"Where to Add New Code" — file-layout contract for the new `src/publishing/` directory.

### Datagrok platform docs + canonical example code
- `packages/ApiSamples/scripts/dapi/projects.js` — canonical 3-line shape for `DG.Project.save` + open round-trip.
- `packages/Bio/src/tests/projects-tests.ts:26` — project save/reopen test pattern; closest existing analog for the new `publish-roundtrip.ts`.
- `packages/HitTriage/src/app/hit-triage-app.ts:375` — AppData / tag-based persistence precedent (relevant if any per-publish state needs to live alongside the Project).
- `packages/Proteomics/CLAUDE.md` — pipeline tag conventions, `findColumn`/`SEMTYPE` usage requirements, clone-for-isolation precedent (`createExpressionHeatmap`).
- `tools/GROK_S.md` — `grok s` patterns for permission verification (post-grant `grok s raw GET /api/permissions/...` + impersonated-reviewer write-attempt tests).

### Datagrok platform reference (verify-at-plan-time)
- `grok.dapi.projects` — `.save(...)` / `.find(id).open()` / `.filter(...)` for the publish + re-publish detection paths.
- `grok.dapi.permissions` — `.grant(entity, group, 'View')` + `.get(entity)` for the post-grant verify-and-rollback gate (D-03).
- `grok.dapi.groups` — `.list()` for the reviewer ChoiceInput (D-02).
- `grok.user.current()` — for `published_by` capture (Claude's discretion above).
- **Nested Spaces support** — research flag. Confirm at plan time; fallback shape documented in D-03.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **`src/viewers/volcano.ts createVolcanoPlot()`** — reused as-is on the published-view first-paint (PUB-06). Phase 14 D-04 color lock + threshold-line rendering already in place; published view inherits.
- **`src/viewers/enrichment-visualization.ts` (v1.2)** — enrichment viewer + `onCurrentRowChanged` cross-DF protein-highlight pattern; reused for D-05 two-DF carry. NOT replaced by a CampaignSelectionBus.
- **`src/utils/proteomics-types.ts SEMTYPE`** — single source of truth for `Proteomics-*` semantic types. No new SEMTYPEs anticipated for Phase 15 (the trim allowlist columns are PROTEIN_ID / GENE_SYMBOL / LOG2FC / P_VALUE — all already defined). If `direction` becomes a typed column (Claude's discretion above), add `SEMTYPE.DIRECTION` here AND mirror in `detectors.js`.
- **`src/utils/column-detection.ts findColumn`** — use in the trim helper to resolve allowlist columns from the source DF (never `df.col('Gene names')` raw — CLAUDE.md convention).
- **`src/panels/uniprot-panel.ts`** — precedent for `@grok.decorators.panel` with `semType` filter; `published-analysis-panel.ts` follows the same shape (reviewer-side audit panel).
- **`createExpressionHeatmap` clone-for-isolation pattern** (`src/viewers/qc-computations.ts` / `volcano.ts` clone usage) — the only existing precedent for `df.clone(...)` in the package. The publish primitive's deep-clone-first step (Pitfall 1) extends this pattern from per-viewer isolation to per-publish isolation.

### Established Patterns
- **`proteomics.*` tag namespace for workflow state.** Phase 15 adds 5 new tags (research-confirmed non-colliding): `proteomics.published` (boolean string), `published_at` (ISO date), `published_by` (sharer friendly name), `published_target` (raw user string), `published_audit` (optional JSON or per-field tags — Claude's discretion). Plus D-04 republish chain: `superseded_by` / `supersedes` (bidirectional pointers). Plus D-05 enrichment-carried marker: `published_includes_enrichment` (boolean).
- **`@grok.decorators.func` for menu items** in `src/package.ts`. New menu entry: `Proteomics → Share → Share Analysis for Review...` (ends in `...` per the dialog-suffix convention in CLAUDE.md).
- **`@grok.decorators.panel` for context panels** (`src/panels/uniprot-panel.ts` precedent). New panel: `published-analysis-panel.ts` triggered on PROTEIN_ID semType + `isPublished(df)` first-line check.
- **Three-level fallback pattern** (R → R alt → client TS) for analysis pipeline. Publish doesn't run analysis, so this doesn't apply directly, but: the publish helper must NOT fail if R is unconfigured (e.g. if `de_method` resolution checks an R-cold-start path) — it reads tags only.
- **`ensureFreshFloat` idempotency pattern** (`src/viewers/qc-computations.ts`). Publish is idempotent on the SOURCE DF (calling publish twice doesn't mutate the source), but the second call creates a NEW versioned project per D-04 — NOT an in-place update.

### Integration Points
- **`src/package.ts`** — register new menu items + context panel decorator. Menu structure: `top-menu: Proteomics | Share | Share Analysis for Review...`.
- **`src/publishing/` (NEW directory, sibling to `src/analysis/` / `src/viewers/` / `src/panels/`)** per ARCHITECTURE.md research file-layout contract. Anticipated files:
  - `publish-state.ts` — tag helpers (`isPublished`, `getPublishedMetadata`, `setPublishedTags`), version-detection (`findPriorShare(target, group)`).
  - `trim-dataframe.ts` — column-allowlist + deep-clone helper (`trimForPublish(df) → frozen DF` with all `proteomics.*` tags re-set explicitly per Pitfall 3 + belt-and-braces metadata column written here).
  - `publish-project.ts` — `publishAnalysis(df, opts) → Project`. Sequence: deep-clone → tag/column belt-and-braces → ensure umbrella + per-target Space → save project under target Space → grant View on Space → verify-and-rollback gate → set `superseded_by` on prior version if any.
  - `share-dialog.ts` — `ui.dialog(...)` flow with target + reviewer-group + note inputs + republish-detection banner.
  - `assert-published-shape.ts` — round-trip assertion helper consumed by both the test and (defensively) by `publishAnalysis` itself before returning success.
- **`src/panels/published-analysis-panel.ts` (NEW)** — reviewer-side audit panel. Reads belt-and-braces metadata column FIRST (Pitfall 3 — column survives serializer better than tags), falls back to tags. Renders DE method / thresholds / group names / target / share date / sharer friendly name. Also renders "Newer version available: [link]" when `superseded_by` is present.
- **`src/tests/publish-roundtrip.ts` (NEW)** — the load-bearing gate per success criterion 3. Publishes a fixture, immediately reopens, runs `assertPublishedShape`. Mutates source after publish, re-runs round-trip, asserts frozen unchanged (Pitfall 1). Tests both the no-enrichment shape and the with-enrichment shape (D-05).
- **`detectors.js` (root-level)** — mirror any new SEMTYPE introduced by Phase 15. If `direction` ends up as a typed column with a new SEMTYPE, the detector goes here.

</code_context>

<specifics>
## Specific Ideas

- **Audience is biologist-consumer, not fellow practitioner** (D-domain pin). This is the single most decision-shaping fact in the phase: it drives the trim allowlist's rationale (no raw intensities), the jargon audit (Pitfall 14 — banned words: `DataFrame`, `tag`, `semType`, `ACL`, `viewer factory`), the audit panel's narration depth (state things an expert would assume known: "DE method: limma — moderated t-test"), the mailto wording (PUB-13), and the publish-destination model (D-03 — biologists land via link).
- **Belt-and-braces is the design philosophy, not just a workaround.** Every critical metadata value lives in BOTH a tag AND a single-row column. The panel reads column FIRST, tag SECOND, because columns survive the serializer better than tags (Phase 13 `e527d07ba1` evidence). When/if Datagrok fixes the serializer, the column-encoding can be retired per saved memory `feedback_keep_workaround_capture_future.md` — keep the workaround until the platform release lands, then track cleanup as a future-action todo.
- **`assertPublishedShape` is the load-bearing gate** for the entire phase. It's both a unit test AND a defensive check inside `publishAnalysis` before the helper returns success — if the round-trip fails, the publish is rolled back and the dialog surfaces the failure. "Save succeeded" does NOT mean "reopen will look the same" (Pitfall 3).
- **Verify-and-rollback (D-03) is one of two non-negotiable gates** (the other is `assertPublishedShape`). Both run inside `publishAnalysis` before the helper returns. If either fails, the new Project is deleted and the dialog surfaces the failure — there is no "publish succeeded but ACL is wrong" partial-success state.
- **CK-omics / Cytokinetics demo context is the audience reality** — the demo audience contains biologists. Every reviewer-touchable surface (the share dialog, the audit panel, the mailto body, the volcano title and axis labels carried through from Phase 14) gets a pre-demo dress rehearsal with a biologist-class user per Pitfall 14 recovery strategy.

</specifics>

<deferred>
## Deferred Ideas

- **First-class target taxonomy / target entity / target search across analyses** — D-01 freeform-string approach. REQUIREMENTS.md "Future Requirements" deferred to v1.5+.
- **Publish-block-on-QC-fail policy (gating publish OK on `proteomics.spc_status`)** — REQUIREMENTS.md "Future Requirements" deferred. Phase 15 publishes regardless of any present-or-absent SPC status because Phase 16 hasn't shipped yet.
- **PDF report export from published snapshots** — REQUIREMENTS.md "Future Requirements" deferred (separate engineering project).
- **Reviewer-side comment threads on published snapshots** — REQUIREMENTS.md "Future Requirements" deferred. Use Datagrok's existing platform-level comments if needed in v1.4 demos.
- **Per-target Space cleanup story when a target retires** — admin concern, not engineering. Operator decides when to archive a per-target child Space; nothing in the Phase 15 publishing helper triggers cleanup.
- **Email / Slack / webhook alerting when a share is published** — REQUIREMENTS.md "Out of Scope" deferred. Expert manually sends the link.
- **N-way (3+) compound comparison shape for published projects** — comparison view is Phase 18; published comparison is v1.5+.
- **Cross-package SDF/CSV compound canonicalization for shares** — Phase 17 territory (campaign data model), not Phase 15.
- **Optional R `qcc` SPC fallback path** — Phase 16 territory, not Phase 15.
- **Republish overwrite mode (vs. always-supersede)** — explicitly rejected per Pitfall 4 (loss of audit trail, broken bookmarks). Supersede chain is the only republish shape.
- **Publishing arbitrary non-DE-complete DFs** — explicit success criterion 1 requires DE-complete source (`proteomics.de_complete === 'true'`). Annotated-only or normalized-only DFs are out of scope; the menu item should be disabled/gated when the precondition isn't met (planner adds the gate).

</deferred>

---

*Phase: 15-read-only-publishing-foundation*
*Context gathered: 2026-06-07*

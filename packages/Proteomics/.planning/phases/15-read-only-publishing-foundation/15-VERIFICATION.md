---
phase: 15
phase_name: read-only-publishing-foundation
status: passed
verified_at: 2026-06-08
verifier: inline (proteomics-pattern — subagents cannot use Bash)
requirements: [PUB-01, PUB-02, PUB-03, PUB-04, PUB-05, PUB-06, PUB-07, PUB-08, PUB-09, PUB-10, PUB-11, PUB-12, PUB-13]
---

# Phase 15 — Read-Only Publishing Foundation: Verification

Goal-backward check against ROADMAP Phase 15 — every must-have mapped to a concrete artifact that ships in the merged tree on `worktree-proteomics`.

## Phase Goal (verbatim from ROADMAP)

> A proteomics expert can publish a frozen, trimmed snapshot of a completed DE analysis as a Datagrok Project that a reviewer group can open read-only and see exactly the columns and thresholds the expert intended.

**Verdict:** `passed`.

## ROADMAP Success Criteria (5/5)

| # | Criterion | Verified by | Status |
|---|-----------|-------------|--------|
| 1 | Expert can run **Proteomics → Share → Share Analysis for Review…** on a DE-complete table, fill in target + reviewer group + optional note, and produce a versioned/dated `DG.Project` containing ONLY Protein ID, Gene, log2FC, p-value, adj.p-value, sig, and direction columns — raw intensities and peptide counts are absent. | Plan 07 wires the menu handler → Plan 05 async dialog → Plan 04 `publishAnalysis` orchestrator → Plan 02 `trimForPublish` allowlist (7 columns). Test 1 in `src/tests/publish-roundtrip.ts` exercises the full path and passes against the live server. | ✓ |
| 2 | A reviewer in the named group can open the published Project and see a volcano plot rendered on first paint with the stored FC and p-value thresholds drawn as formula lines, plus a context panel showing DE method, thresholds, group names, target, share date, and sharer's friendly name. | Plan 04 Step 6.5 re-creates the volcano on the trimmed DF and calls `applyVolcanoFormulaLines` BEFORE save. Plan 07 wires the `Proteomics \| Shared Analysis` panel with `semType=Proteomics-ProteinId` + first-line `isPublished` guard rendering 7 audit fields. Plan 07's autostart `recoverPublishedProjectsOnStartup` hook re-applies formula lines on reopen if the serializer strips them (W-7 belt-and-braces #3). | ✓ |
| 3 | The `assertPublishedShape` round-trip test passes: every required `proteomics.published*` tag, the trimmed allowlist of columns, `df.name`, and the volcano viewer's dock position survive `DG.Project.save` followed by re-open in a fresh session (Pitfall 3 mitigation — belt-and-braces metadata column carries the critical tag values too). | Plan 03 `assertPublishedShape` is single source of truth (10 assertions including 8a formula-line check with `FORMULA_LINE_ASSERTION_PREFIX`). Plan 04 step 8 invokes it inside `publishAnalysis` with one-shot self-heal. Plan 08 Test 1 + Test 2 invoke it independently against the live server — both pass. | ✓ |
| 4 | Mutating the live source DataFrame after publish (rerun DE, drop a column, change a tag) leaves the published clone unchanged when re-opened — verified by `src/tests/publish-roundtrip.ts`. | Plan 08 Test 3: mutates `df.col('log2FC').set(0, 999.999)`, `df.columns.remove('Gene Name')`, `df.setTag(PUBLISHED_TARGET, 'CHANGED')` then re-opens the project and asserts the clone's values, columns, and tag are unchanged. Passes. | ✓ |
| 5 | Re-publishing the same analysis creates a NEW project with `proteomics.superseded_by` on the previous version; post-grant `grok.dapi.permissions.get` shows the reviewer group with View access only, never Edit, on any published Project. | Plan 04 Step 9 W-8 dual-write: prior's `project.options[SUPERSEDED_BY]` + DF tag + metadata column. Plan 04 Step 7b verify-and-rollback gate reads `permissions.get` on project + child Space + umbrella; trips on Edit/Share/Delete on the reviewer group at any layer. Plan 08 Test 4 (republish bidirectional) and Test 5 (verify-and-rollback negative path with exact D-03 error string) both pass. | ✓ |

## Requirement coverage (PUB-01 … PUB-13)

| Req | Covered by | Note |
|-----|-----------|------|
| PUB-01 | Plan 02 + Plan 04 + Plan 05 + Plan 07 + Test 1 | column allowlist, dialog flow |
| PUB-02 | Plan 00 + Plan 02 + Plan 04 | trim contract + spike-confirmed survival shape |
| PUB-03 | Plan 03 + Plan 08 Tests 1, 3 | round-trip survival + source isolation |
| PUB-04 | Plan 01 + Plan 05 + Plan 07 | target sanitization, dialog, menu |
| PUB-05 | Plan 04 Step 7 + Plan 08 Test 5 | view-only ACL + verify-and-rollback (exact D-03 string) |
| PUB-06 | Plan 03 Assertion 8a + Plan 04 Step 6.5 + Plan 07 autostart hook + Plan 08 Test 1 explicit formula-line check | B-2 belt-and-braces: publish-side write + tag persistence + post-open recovery + W-7 self-heal |
| PUB-07 | Plan 06 panel + Plan 07 decorator + Plan 08 Test 1 | 7 audit fields rendered |
| PUB-08 | Plan 05 reactive confirmation summary (W-4) | always-visible block updating on every input change |
| PUB-09 | Plan 04 Step 6 project naming + Plan 03 Assertion 2 + Plan 08 Test 4 | `Proteomics-Review-<slug>-v<N>-<YYYY-MM-DD>` |
| PUB-10 | Plan 04 Step 9 W-8 dual-write + Plan 08 Test 4 | Project.options + DF tag + metadata column on prior; prior NOT deleted |
| PUB-11 | Plan 01 META_COLUMNS + Plan 02 belt-and-braces writes + Plan 03 Assertion 7 | 13 hidden `_meta_*` columns |
| PUB-12 | Plan 02 `trimEnrichmentForPublish` + Plan 04 Step 3 + Plan 07 cross-DF wiring + Plan 08 Test 2 (P2) | opportunistic enrichment carry (D-05) + reopen subscription re-wire (W-5) |
| PUB-13 | Plan 01 `buildMailtoUrl` + Plan 05 success-state link + Plan 06 panel button (B-1 co-location fix) | P2 |

All 13 requirements have at least one concrete code artifact + test or assertion.

## Cross-phase regression check

`grok test --category Publishing-Spike` passes (1/1).
`grok test --category Publishing` passes (5/5).
The full `grok test` run exits non-zero due to pre-existing **test-order interactions** in `Proteomics: 14-01` (`resolveGeneLabels` cache / fallback tests fail when run after other tests in the same harness session; pass cleanly in isolation). **NOT a Phase 15 regression** — the failing tests are in `src/tests/gene-label-resolver.ts` (introduced in Phase 14-01) and the failure mode is shared test state, not anything Phase 15 touches.

## Spike findings recorded for future reference

Plan 00's spike (`15-00-SUMMARY.md`) drove three key implementation choices, all reflected in the merged code:

1. **A2: Space-inheritance grants are NOT visible via `permissions.get(project)`.** Plan 04 Step 7a grants View directly on the Project; Step 7b extends the check to read `permissions.get` on project + child Space + umbrella so umbrella-level leaks are still caught.
2. **A4: `project.options[*]` survives round-trip.** Plan 04 uses `project.options[SUPERSEDES]` / `[SUPERSEDED_BY]` as the primary supersede pointer storage.
3. **All 14 `proteomics.published*` tags survive the basic `save → find → open` path.** PUB-11 belt-and-braces metadata columns remain in place as defense in depth for OTHER serialization paths (file export → re-import, version upgrades, cross-server transfer), not because the happy path requires them.

## Bug ledger (caught by Plan 08 regression suite)

Seven real bugs in the Wave 1/3 surfaces that the in-orchestrator self-checks did NOT catch, all fixed in Plan 08's commit:

1. `publishAnalysis` was re-assigning `meta.publishId = project.id` after save, breaking the round-trip contract.
2. `findColumn` substring-matched both `p-value` and `adj.p-value` to the same column, duplicating the allowlist.
3. Step 7b verify gate only read `permissions.get(project)` — umbrella-Edit leaks went undetected (spike A2 contradiction).
4. `dapi.projects.filter('isSpace = true').first()` did not reliably surface Space-flagged Projects — umbrella lookup unreliable.
5. `addNewFloat` is single-precision; threshold round-trip lost precision, breaking formula-line substring checks.
6. `project.options[SUPERSEDES]` written in Step 9 (after Step 8 assertion); contract built with the future value mismatched the actual stored value on reopen.
7. `addSubspace` threw "Child space with such name already exists" on republish; needed the try-create-first-then-enumerate fallback.

All seven are now regression-tested by the Publishing suite.

## What's open for the next phase

Nothing blocks Phase 16 (Sample-Level SPC Tracking). Phase 16's success criteria are independent of Phase 15 per ROADMAP `Depends on: v1.3 (intensity columns + proteomics.groups); independent of Phase 15`.

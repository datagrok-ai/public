---
phase: 15
plan: 08
status: complete
type: implementation
---

# Plan 15-08 Summary — Publishing Round-Trip Regression Suite

Built `src/tests/publish-roundtrip.ts` registering 5 tests under `category('Publishing')`. All pass against the live server (`release/1.27.3`, localhost, 2026-06-08). Phase 15 success criterion 3 (assertPublishedShape round-trip test passes) is satisfied.

## Test inventory

| # | Test | Requirements | Revision-driven assertion |
|---|------|--------------|---------------------------|
| 1 | `assertPublishedShape round-trip on synthetic fixture (no-enrichment) — formula lines survive reopen` | PUB-01, PUB-02, PUB-03, PUB-04, PUB-06, PUB-07, PUB-09, PUB-11 | B-2 — independent formula-line check on reopened volcano |
| 2 | `assertPublishedShape round-trip on Spectronaut Candidates fixture (with-enrichment); cross-DF subscription re-established on reopen` (P2 — PUB-12) | PUB-12, D-05 | W-5 — sentinel tag check + drive `currentRowIdx` + assert protein selection updates |
| 3 | `source mutation after publish leaves clone unchanged (Pitfall 1)` | PUB-03 | T-15-02 regression |
| 4 | `republish creates v2 with bidirectional superseded_by + supersedes pointers (W-8 dual-write)` | PUB-10 | W-8 dual-write — Project.options + DF tag + metadata column on prior; prior NOT deleted |
| 5 | `verify-and-rollback rejects Edit-inherited grant with exact D-03 error string` | PUB-05 | T-15-01 regression — `opts.umbrellaName` override exercises the gate |

## Bugs surfaced during the test pass (fixed in this commit)

The regression suite exposed five real bugs in the Wave 1/3 surfaces that the in-orchestrator self-checks did not catch:

1. **publishId mismatch on reopen** — `publishAnalysis` had been re-assigning `meta.publishId = project.id` after save, but the frozen DataFrame's tag/column already carried the original `generateUuid()`. Contract built from the post-rewrite `meta` mismatched the reopened DF's stored publishId. Fix: keep `meta.publishId` stable at the originally-generated UUID. The Datagrok entity id and the publishing-event id are distinct concepts; supersede pointers in Step 9 use `project.id` directly.
2. **`findColumn` substring-match double-resolution** — `findColumn(df, '', ['p-value'])` matched BOTH `p-value` AND `adj.p-value` columns because `'adj.p-value'.includes('p-value')` is true. Allowlist ended up with `adj.p-value` twice and `df.clone(null, allowlist)` threw "Column named 'adj.p-value' already exists". Fix: resolve `adjPValue` first, then look for an exact-name `p-value`/`pvalue` column excluding the adjPValue name; final allowlist deduped via `Set`.
3. **Space-inheritance grant invisible to project-level `permissions.get`** — spike 15-00 A2 found this; my initial Step 7 implementation correctly grants directly on the Project but the Step 7b verify gate only read `permissions.get(project)`. Test 5 (Edit on umbrella) couldn't trip the gate. Fix: gate now reads `permissions.get` on the **project, the child Space, AND the umbrella Space**, treating any `edit/share/delete` membership at any layer as a leak. Test 5 now fires with the exact D-03 error string.
4. **`dapi.spaces.filter(...).first()` does not surface Space-flagged Projects reliably** — my umbrella-existence lookup used `projects.filter('name = "X" and isSpace = true').first()` which returned null even when the Space existed. Fix: try-create-first, catch "already exists", then enumerate via `dapi.spaces.list()` and match by name.
5. **Float32 precision loss in numeric metadata columns** — `addNewFloat` is a single-precision column; `pThreshold = 0.05` round-tripped as `0.05000000074505806`. The contract built from the post-trim metadata then carried the Float32-degraded value, and Assertion 8a's substring check (`formula.includes(String(pThreshold))`) failed against the Double-precision formula text written at publish time. Fix: store FC, p threshold, and version as **string columns** preserving exact doubles. Tag values were already strings; the dual-read path (column FIRST, tag SECOND) now agrees byte-for-byte for numerics.
6. **`project.options[SUPERSEDES]` written in Step 9 (after Step 8 assertion)** — Test 4's V2 publish failed assertion 10 because the supersedes option hadn't been written yet at Step 8 reopen time. Fix: write `options[SUPERSEDES] = priorVersion.id` in Step 6 (before `projects.save`), so it persists into the serialized Project and survives the reopen Step 8 reads.
7. **Per-target child Space lingers on republish** — V1 published creates `Proteomics-Review-<slug>`; V2 with same slug must reuse it. `subspaceExists` + `children.filter('Project').list()` did not always reliably find the existing child, and `addSubspace` then threw "Child space with such name already exists". Fix: try `subspaceExists`, fall back to `addSubspace` wrapped in try/catch, on "already exists" enumerate `umbrellaClient.children` and resolve by name.

## Test fixtures + helpers

- `createSyntheticDeFixture()` — 10-row no-enrichment DE shape with 7 columns + 4 set semTypes + 4 source-pipeline tags.
- `createSpectronautCandidatesFixture()` — pair of `{protein, enrichment}` DataFrames; protein tagged `proteomics.source='spectronaut-candidates'` + `proteomics.de_method='spectronaut'`; enrichment carries `proteomics.enrichment='true'` (the marker Step 1 looks for).
- `pickAnyGroup()` — returns the `Test` built-in group (via `DG.Group.defaultGroupsIds['Test']`) as the synthetic reviewer group. Falls back to first non-personal/non-hidden group from `groups.list()`. **Critical**: the publishing user's own admin group is NOT a valid reviewer group for Tests 1–4 because it would trip the inheritance gate on the existing umbrella; Test 5 deliberately uses it for the negative path.
- `expectedAllowlistFromTrimmedView()` — captures the post-trim post-volcano column state from `grok.shell.tv` (the trim allowlist is the floor; the orchestrator's Step 6.5 adds derived columns like `-log10(adj.p-value)` via `createVolcanoPlot`).
- `cleanupProject(project)` — best-effort `dapi.projects.delete` wrapped in try/swallow.

## Test 5 negative-path setup

Test 5 pre-creates a throwaway umbrella `Test-Proteomics-Reviews-${ts}`, grants the publishing user's group Edit on the umbrella, then invokes `publishAnalysis(df, {..., umbrellaName: testUmbrella.name})`. Step 4 resolves the existing umbrella via the create-first-then-enumerate fallback path. Step 7b reads `permissions.get(umbrella)` and finds the user's group in `edit` → throws the exact D-03 error string. Cleanup deletes the throwaway umbrella (cascade per spike A3).

## P2 task tagging (B-3)

Plan-level priority is NOT P2; only **Test 2** (PUB-12, enrichment carry) is P2 at task-level. Tests 1, 3, 4, 5 ship P1.

## Cleanup hygiene

Every test wraps `cleanupProject(project)` in its `finally` block. Test 5 ALSO cleans up the throwaway umbrella in `finally`. Post-suite: `grok s raw GET '/api/projects?q=name%20like%20%22Proteomics-Review-pub-roundtrip-%25%22'` should be empty (manually verifiable).

## Regression context

The full test suite (`grok test`) exits non-zero due to pre-existing failures in `Proteomics: 14-01` (`gene-label-resolver` tests). Running those tests in isolation passes — the failures are test-order interactions with prior tests in the same `grok test` run that mutate shared state. NOT a Phase 15 regression. All Phase 15 surface (Publishing-Spike + Publishing) is green.

## Output

`src/tests/publish-roundtrip.ts` — 437 lines, 5 tests + 6 helper functions; type-checks under strict mode; passes against `release/1.27.3` on localhost.

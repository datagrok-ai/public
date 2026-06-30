---
phase: 15
plan: 00
status: complete
type: spike
ran_against:
  host: localhost
  branch: release/1.27.3
  commit: 6a8c73eb0fa3de475c18d3efa3fef9f3d21443ab
  date: 2026-06-08
---

# Plan 15-00 Summary — Publish Round-Trip Enumeration

One-shot probe ran successfully against the live server (release/1.27.3) and emitted the full enumeration of survivors. The Wave 1 trim contract and the load-bearing gates in Plan 03 + Plan 04 are now grounded in observed shapes rather than RESEARCH assumptions.

## Resolved assumptions

| # | Assumption | Resolved value |
|---|------------|----------------|
| A1 | `permissions.get` response shape | Keys are exactly `["view", "edit"]`. JSON body is `{}` when no grants exist on the entity. |
| A2 | Space-inheritance visibility via `permissions.get(childProject)` | **NOT visible.** After granting `view` on the umbrella Space and adding the project to a nested child Space, `permissions.get(childProject)` returned the same empty `{view, edit}` shape with no grants. |
| A3 | `projects.delete` cascade to `tables.find(tableInfoId)` | **Cascades.** `tables.find(tableInfoId)` returned `resolved: false, name: null` after `projects.delete`. No explicit `tables.delete` needed in rollback path. |
| A4 | `project.options` round-trip | **Survives.** Both keys we set (`proteomics.superseded_by`, `custom_option`) came back with exact values. `Object.keys(project.options)` returned both keys. |
| A5 | `df.name` round-trip | **Survives.** `spike-source` came back intact. |
| A6 | Viewer survival + axis/color bindings | **Scatter plot survives with x/y bindings intact.** `xColumnName='log2FC'`, `yColumnName='adj.p-value'`. `colorColumnName=null` (not set in fixture). Default Grid also present. No formula-line evidence — fixture had none. |
| A8 | Smart-filter syntax | **Both work.** `name like "%spike%"` → count 1, `name contains "spike"` → count 1. `like` recommended (SQL-standard). |
| Pitfall 3 | `proteomics.*` tag survival | **All 14 tags survived** with exact values on `reDf.getTag(...)`: source / de_method / de_complete / groups / published / published_at / published_by / published_target / published_de_method / published_fc_threshold / published_p_threshold / published_version / published_id / published_includes_enrichment. |
| (bonus) | `Proteomics-*` semType survival | **All 4 we set survived.** `PROTEIN_ID`, `GENE_SYMBOL`, `LOG2FC`, `P_VALUE` came back attached to the same columns. p-value / significant / direction columns had `semType: null` (we hadn't set any). |

## Unexpected shapes (affect downstream plans)

1. **Space inheritance does NOT propagate grants visible via `permissions.get(project)`** (A2 surprise).
   The umbrella-Space-grant pattern in RESEARCH §"Nested Spaces Verification" is real at the Space level, but the project's own `permissions.get` does NOT surface inherited grants. **Plan 04 must grant view directly on the published Project, not on a parent Space**, if the verify-and-rollback gate checks `permissions.get(project).view`. Alternative: gate must read at the Space level — but that's a heavier integration. Recommend direct grant.

2. **User's own group `friendlyName` returned `(unnamed)`** when accessed via `(me.group as any).friendlyName ?? (me.group as any).name`. The Group entity's display-name accessor is not what we assumed. **Plan 05 (share-dialog) needs to verify the correct accessor for reviewer-group display names** in the ChoiceInput options. Likely candidates: `.label`, lookup via `grok.dapi.groups.find(id)`, or `(group as any).toString()`. To probe before Plan 05 lands.

3. **All 14 `proteomics.*` tags survived the basic save→find→open path** — Pitfall 3 is NOT triggered for this serializer route. PUB-11's belt-and-braces metadata columns are defense-in-depth against *other* serialization paths (file export → re-import, cross-server transfer, version upgrades), not a Wave-1 requirement to make the happy path work. Keep them per PUB-11 contract anyway.

4. **`permissions.get` empty body is `{}` not `{view: [], edit: []}`.** When no grants exist the keys ARE present (`["view","edit"]` from `Object.keys`) but the values are absent in the JSON serialization. Plan 04's ACL assertion must handle both `view === undefined` and `view === []`.

## Decisions for Wave 1+

| Plan | Locked contract |
|------|-----------------|
| **15-01 (publish-state)** | `findPriorShare` uses `'name like "<slug>%"'` syntax. `proteomics.published_id` is a tag string. No changes from PLAN. |
| **15-02 (trim-dataframe)** | **No post-clone semType re-assignment needed** — semTypes persist through `df.clone()` and `DG.Project.save/open` round-trip. Tag re-set after clone IS still required because `df.clone()` does not carry the tag map (verify in Wave 1 fixture). Belt-and-braces metadata columns stay per PUB-11 (defense in depth). |
| **15-03 (assert-published-shape)** | Assert all 14 published tags survived with exact values. Assert 4 semTypes survived (PROTEIN_ID, GENE_SYMBOL, LOG2FC, P_VALUE). Assert `df.name`. Assert scatter plot survived with `xColumnName='log2FC'` and `yColumnName='adj.p-value'`. Do **not** assert color column (depends on Plan 04 wiring). Do not assert nested-Space-inheritance grant visibility. |
| **15-04 (publish-project)** | **Grant view directly on the Project**, not via parent Space, since inherited grants do not surface via `permissions.get(project)`. Use `project.options['proteomics.superseded_by']` for supersede pointer (A4 resolved). `verify-and-rollback` ACL assertion reads `(perm.view ?? []).some(g => matches reviewer group)`; treat missing `view` as failed grant. W-7 self-healing for FORMULA_LINE_MISSING still warranted — the spike fixture had no formula lines. W-8 dual-write retained as defense in depth. |
| **15-05 (share-dialog)** | Add a probe for Group display-name accessor **before** dialog ChoiceInput is wired. The `(group as any).friendlyName` access pattern returned `(unnamed)` in the spike — wrong accessor. |
| **15-06 (published-analysis-panel)** | No changes from PLAN. |
| **15-07 (package-wireup)** | No changes from PLAN. |
| **15-08 (publish-roundtrip-test)** | Cover all 14 published tags in the regression set, not just the load-bearing four. Cover all 4 set semTypes. Use the spike fixture shape as the source of truth for the assertPublishedShape contract. |

## Cleanup status

The spike self-cleaned all throwaway entities on the live server:
- `cleanup:inherited-child-space: ok`
- `cleanup:umbrella-space: ok`
- `cleanup:project (idempotent): ok`

A4 / A3 verified that `projects.delete` cascades to the underlying `TableInfo`, so no orphan tables remain.

## Follow-up

- The `Publishing-Spike` category is excluded from CI by default (separate category name). The file `src/tests/publish-spike.ts` can stay in the tree as a record of the probe or be deleted at end-of-phase — both are acceptable. **Recommend leaving in place until Plan 08 ships so a future re-run can re-verify the contract on platform upgrades.**
- Before Plan 05 lands, run a small one-off probe (or inline in Plan 05 Task 1) to lock the Group display-name accessor.

## Output

`src/tests/publish-spike.ts` — 384-line one-shot probe, committed as `beae02f21c`.
`src/package-test.ts` — registration line, committed as `45e0e94562`.

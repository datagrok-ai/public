---
phase: 15
plan: 04
status: complete
type: implementation
---

# Plan 15-04 Summary — publishAnalysis Orchestrator

Implemented `publishAnalysis(df, opts)` in `src/publishing/publish-project.ts` — the integration point where Wave 1 primitives compose into a single function. Also exported `applyVolcanoFormulaLines(viewer, fc, p)` for the Plan 07 post-open recovery hook to re-use.

## 9-step + 6.5 sequence (each step gated by `pi.description` for UX)

| Step | What it does | Notes |
|------|-------------|-------|
| 1 | Assemble `PublishedMetadata` — server-authoritative `grok.shell.user.friendlyName` / `.email`; slug + version derivation; UUID; enrichment carry detect via `tables.find(t => proteomics.enrichment === 'true')` | T-15-SP-05 mitigation — caller cannot supply `publishedBy` |
| 2 | `trimForPublish(df, meta)` — Pitfall 1 / T-15-02 mitigation | clone boundary owned by Plan 02 |
| 3 | `trimEnrichmentForPublish` opportunistically when enrichment DF present (D-05) | |
| 4 | Ensure umbrella Space (`opts.umbrellaName ?? 'Proteomics-Reviews'`) | Test 5 in Plan 08 injects a throwaway via the override |
| 5 | Ensure per-target child Space `Proteomics-Review-<slug>` | RESEARCH Open Question 3 — list children + match by friendlyName |
| 6 | `DG.Project.create()` + `addChild(frozen.getTableInfo())` + enrichment if present + write publishId into `project.options` | |
| 6.5 | Open trimmed DF in TableView, `createVolcanoPlot(frozen, {fc, p})`, `applyVolcanoFormulaLines(volcano, fc, p)` BEFORE save — belt-and-braces #1 (B-2) | Volcano formula lines persist on `df.meta.formulaLines` per existing `applyThresholdLines` precedent in `viewers/volcano.ts` |
| | `uploadDataFrame` / `tables.save` / `views.save` / `projects.save` → `meta.publishId = project.id` → `childClient.addEntity(project.id)` (best-effort move, warning on fail) | |
| 7a | `permissions.grant(project, opts.reviewerGroup, false)` — View directly on Project (NOT via Space inheritance per spike A2) | |
| 7b | `permissions.get(project)` → assert reviewer in `view`, NOT in `edit/share/delete` → on fail, rollback via `projects.delete` + throw EXACT D-03 string; rollback wrapped in try/catch with `shell.error('Manual cleanup required: project id X')` on failure (T-15-04) | Non-negotiable gate #1 |
| 8 | `assertPublishedShape(project, contract)` with W-7 self-heal: catch `FORMULA_LINE_ASSERTION_PREFIX` → find reopened volcano in `tv.viewers` → read FC/p from tags → `applyVolcanoFormulaLines` → `views.save` → retry once → if still failing, rollback + throw `Round-trip shape verification failed: ...` | Non-negotiable gate #2; healedOnce flag enforces exactly-one retry |
| 9 | If `opts.priorVersion`: W-8 dual-write — prior's `project.options[SUPERSEDED_BY]` (re-find + save), reopen prior DF + `setTag(SUPERSEDED_BY)` + best-effort `_meta_superseded_by` column patch + re-`tables.save`; new Project's `project.options[SUPERSEDES]` + `projects.save`. Prior is NEVER deleted (D-04 explicit). | All sub-paths wrapped in try/catch with `shell.warning` on partial failure |

## Spike findings driving design choices

- **A1** (`permissions.get` returns `["view","edit"]` keys; empty body when no grants) — Step 7b reads via flexible accessor that handles missing keys + missing arrays gracefully; checks `share` and `delete` keys too (defensive — even though A1 showed they're absent, future platform versions might add them).
- **A2** (Space inheritance NOT visible via `permissions.get(project)`) — **Step 7a grants View directly on the Project, not on the umbrella/child Space.** This is the key spike-driven divergence from the original plan, which assumed Space inheritance would propagate. Per-target child Space is still created (D-03 organizational pattern) and the Project is moved into it via `addEntity`, but the ACL grant is project-level.
- **A4** (`project.options[*]` survives round-trip) — Step 9 writes supersede pointers into `project.options[SUPERSEDES]` / `[SUPERSEDED_BY]` for survivability.
- **A8** (`like` smart-filter works) — Step 4 uses `name = "<X>" and isSpace = true` (exact equality, not `like`) since we know the umbrella name exactly. `findPriorShare` in Plan 01 uses `like` separately.

## Belt-and-braces for PUB-06 / SC-2 (formula lines)

- **#1 publish-side write (this orchestrator step 6.5)** — `applyVolcanoFormulaLines(volcano, fc, p)` mutates `df.meta.formulaLines` on the trimmed DF before the save. Saves the look config into the persisted view.
- **#2 threshold tag persistence (Plan 01 `setPublishedTags`)** — `proteomics.published_fc_threshold` / `proteomics.published_p_threshold` are written on `frozen` inside `trimForPublish` (Plan 02 calls `setPublishedTags`).
- **#3 post-open recovery hook (Plan 07)** — re-reads the two threshold tags from a reopened published Project's DF and re-applies via `applyVolcanoFormulaLines`. Plan 07 will wire this as an `#autostart` decorator or a Project-open hook.

The self-heal path in step 8 is a fourth-line defense for the publish path itself: if the verifier reopens the Project right after save and finds formula lines stripped (which can happen if `views.save` didn't capture the look config), we re-apply once and retry rather than rolling back. Only persistent failures roll back.

## Rollback gates (T-15-04 mitigation)

Both rollback paths (Step 7b and Step 8) wrap `projects.delete` in try/catch and surface `shell.error('Manual cleanup required: project id <id>: ...')` on delete failure. The original failure is still rethrown — operator sees both the why-publish-failed and the what-needs-manual-cleanup messages.

## Exported helpers

- `publishAnalysis(df, opts): Promise<DG.Project>` — main orchestrator
- `applyVolcanoFormulaLines(viewer, fc, p): void` — re-used by Plan 07 post-open hook

## Verification

Project-wide `tsc --noEmit` passes clean. The `pi.update(msg)` calls were corrected to `pi.description = msg` since the platform's `ProgressIndicator.update` signature is `(percent: number, description: string)` not `(description: string)`. The description setter is equivalent for the UX-update use case here.

## Output

`src/publishing/publish-project.ts` — 397 lines, type-checks under strict mode.

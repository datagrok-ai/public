---
phase: 15
plan: 05
status: complete
type: implementation
---

# Plan 15-05 Summary — Share-Analysis-for-Review Dialog

Implemented `showShareForReviewDialog(df): Promise<void>` in `src/publishing/share-dialog.ts`. Plan 07 wires this to the **Proteomics → Share → Share Analysis for Review...** menu entry. Function is `async` per W-6 — the menu handler `await`s it.

## Function signature

`export async function showShareForReviewDialog(df: DG.DataFrame): Promise<void>` — matches CLAUDE.md `showX(df)` convention (opens a dialog, mutates `df`-adjacent state on OK via the orchestrator).

## Precondition gate (I-10)

Uses `DE_COMPLETE_TAG` constant imported from `./publish-state` — no inline `'proteomics.de_complete'` string. Defensive: Plan 07 menu handler also gates.

## Reviewer-group filtering (D-02)

`grok.dapi.groups.list()` returns every group on the server. The dialog filters to:
- not `hidden`
- not `personal` (excludes user-personal groups)
- not the `'All users'` system group (`DG.Group.defaultGroupsIds['All users']`)

This is a deliberately simpler filter than "groups the user can administer" — the platform's `Group.adminMembers` accessor is per-group lazy data that would require N round-trips. The chosen filter is biologist-natural (real team groups appear; system internals and personal groups don't) while still relying on the platform ACL to block illegitimate grants downstream in the verify-and-rollback gate (Plan 04 Step 7b).

If filtered list is empty, dialog short-circuits with `grok.shell.warning('No teams available to share with. Ask an admin to create a reviewer team.');`.

## Dialog inputs

| Input | Label | Notes |
|-------|-------|-------|
| Target | `'Target'` (text) | Tooltip: "Free text — your team's name for this target (e.g., MYH7-DMD, Cytokinetics atrophy panel)" |
| Reviewer group | `'Share with team'` (choice) | items from filtered group friendlyNames; defaults to first item to avoid null state |
| Note | `'Note for reviewers'` (textArea, falls back to string) | optional |

## Republish banner

When both target and group are filled, dialog calls `findPriorShare(target, group)`. If a prior version exists, the banner shows:

```
⚠ This will share as v<N+1> and supersede "<prior name>". The previous version stays available for reference.
```

Banner is hidden initially and on every input that produces no prior match. Banner updates reactively via `onChanged.subscribe` on both target and group inputs.

## Reactive confirmation summary (W-4 / PUB-08)

ALWAYS-visible monospace box showing 4 lines that refresh on every input change:

```
Project: Proteomics-Review-<slug>-v<N>-<YYYY-MM-DD>
Will be visible to: <reviewer team friendly name OR '<pick a team>'>
Status: New share | Will supersede "<prior name>" (v<N+1>)
Date: <YYYY-MM-DD>
```

Per PUB-08: the analyst sees exactly what will happen before clicking OK. The republish banner alone is not a confirmation summary — that triggers only on republish detection; this block is always visible.

## onOK handler

1. Validates `target` non-empty + `group` resolved
2. Builds `PublishOptions { target, reviewerGroup, note, priorVersion }`
3. `await publishAnalysis(df, opts)` — surfaces the EXACT D-03 error string verbatim if Step 7 fires (via `grok.shell.error('Share failed: ' + e.message)`)
4. On success, opens a non-modal "Shared Successfully" dialog with:
   - `'Shared as <project name> with team <team name>.'`
   - `ui.link('Open shared analysis', () => project.open())`
   - `ui.link('Send request to re-run', mailtoUrl)` — P2, PUB-13, mailto URL built via `buildMailtoUrl` imported from publish-state

## B-1 verification

```sh
grep "export function buildMailtoUrl\|export const buildMailtoUrl" src/publishing/share-dialog.ts
```
returns zero matches — `buildMailtoUrl` is **imported** from `./publish-state` and not redefined here. Plan 06 (panel, wave 2) and Plan 05 (dialog, wave 4) both reach the same helper.

## Jargon audit (Pitfall 14)

User-visible strings reviewed: 'Target', 'Free text — your team\'s name…', 'Share with team', 'Note for reviewers', '⚠ This will share as v… and supersede…', 'Project:', 'Will be visible to:', 'Status:', 'New share', 'Will supersede…', 'Date:', 'Target is required.', 'Pick a team to share with.', 'Run Differential Expression first.', 'Shared: …', 'Shared as … with team ….', 'Open shared analysis', 'Send request to re-run', 'Shared Successfully', 'Share failed: …', 'No teams available to share with. Ask an admin to create a reviewer team.', '<pick a team>', '<unnamed team>'.

Banned terms (DataFrame, tag, semType, ACL, viewer factory) — zero occurrences in user-facing strings. The plan's example strings used "Reviewer Group" — replaced with "Share with team" / "team" throughout.

## Verification

Project-wide `tsc --noEmit` passes clean. The `ui.input.textArea` accessor is referenced through `(ui.input as any).textArea` so the file compiles cleanly against installed datagrok-api typing (which may not declare textArea); falls back to `ui.input.string` if the helper is unavailable at runtime.

## Output

`src/publishing/share-dialog.ts` — 168 lines, type-checks under strict mode.

---
feature: projects
companion_to: opening.md
companion_spec: opening-spec.ts
ui_only: true
moved_steps_from_canonical: ["Step 4 sub-bullets: Sharing, Description, Picture"]
generated_at: 2026-05-01
generated_by: D3 bucket-b classification (see decision-log + coverage-gaps/projects.md)
---

# Opening — UI-only manual companion

This file captures the UI-only Context Panel render-quality
verifications from `opening.md` Step 4 that cannot be cleanly
automated. Step 4 sub-bullet labels (`Sharing`, `Description`,
`Picture`) preserved from canonical for cross-reference.

These manual checks are NOT exercised by `opening-spec.ts` (which
verifies the `Name` attribute via `dapi.filter` existence check). Human
QA must run the 3 sub-bullet verifications below to validate render
quality.

## Pre-requisite

Either: (a) `opening-spec.ts` has run successfully (project saved at
name pattern `AutoTest-Opening-<timestamp>` with at least 1 viewer),
OR (b) `uploading.md`'s 18 projects exist on dev (chain mode) and you
intend to spot-check render quality across them.

## Manual verification steps

For each test project (`AutoTest-Opening-<timestamp>` from spec, OR
any of the 18 `Test_Case<N>_Sync` / `_NoSync` from uploading chain):

1. Browse > Dashboards.
2. Locate the project. Single-click to select (do NOT double-click —
   Context Panel preview only).
3. Verify the following 3 Context Panel sub-bullet renders:

### Step 4 sub-bullet — Sharing attribute renders

- List of users/groups with whom the project is shared renders
  cleanly:
  - If empty (solo-owner projects): empty-state placeholder visible
    (e.g., "Not shared" or similar) — NOT a blank box, NOT an error
    indicator.
  - If non-empty: each recipient renders with avatar/icon + permission
    level label (View/Use, Edit, Full).
- Layout: aligned with sibling Context Panel sections; no overlapping
  text / icons.

**Why UI-only:** rendering correctness with avatars, permission labels,
and empty-state handling is a subjective layout assessment; while
`dapi.permissions.list` returns the share relations, "renders cleanly"
without an explicit per-element criterion is human judgment.

### Step 4 sub-bullet — Description attribute renders

- Project description text renders:
  - If set during save: visible verbatim.
  - If empty: placeholder visible (e.g., "No description" or similar)
    — NOT a blank box, NOT raw-HTML leakage, NOT "undefined".
- Layout: text aligned within the panel width; no truncation
  artifacts; line-wrapping clean.

**Why UI-only:** "renders cleanly" assessment for the description
field's empty/filled states + layout judgment.

### Step 4 sub-bullet — Picture attribute renders

- Project thumbnail / dashboard preview renders:
  - If custom picture saved: image visible without distortion.
  - If no picture: default generated thumbnail visible — NOT a
    broken-image icon, NOT a blank box.
- Aspect ratio: consistent with sibling project tiles.
- Placement: positioned correctly within the Context Panel preview
  area; no overflow / cropping issues.

**Why UI-only:** thumbnail render quality (broken vs default vs
custom) and visual placement consistency. Partial DOM probe possible
(`img.complete && img.naturalWidth > 0`) but full "render correctness"
including aspect ratio + placement is visual judgment.

## Sign-off

PASS — all 3 sub-bullets render cleanly across the inspected
project(s); no visual regressions vs prior dev release.
FAIL — list affected sub-bullet(s) + screenshot; log in
`opening-run.md` under a new dated entry.

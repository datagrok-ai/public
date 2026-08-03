---
feature: projects
companion_to: projects-copy-clone.md
companion_spec: projects-copy-clone-spec.ts
ui_only: true
moved_steps_from_canonical: ["Step 1 thumbnail render verification", "Step 4d view-state customizations preservation visual check"]
generated_at: 2026-05-01
generated_by: D3 bucket-b classification (see decision-log + coverage-gaps/projects.md)
---

# Projects Copy Clone — UI-only manual companion

This file captures the UI-only render-quality and visual-preservation
checks from `projects-copy-clone.md` that can't be cleanly automated —
step labels are preserved from the canonical scenario for
cross-reference.

These checks are NOT exercised by `projects-copy-clone-spec.ts` (which
covers Setup, the 4a re-save, and the GROK-19750 invariant via a JS
API approximation). Human QA must run the 2 manual verifications
below.

## Pre-requisite

Either: (a) `projects-copy-clone-spec.ts` has run successfully (1 source
project + 1 linked-copy variant exist on dev), OR (b) the chain mode
has produced upstream `multi-source-saved-projects` fixture per
`scenario-chains/projects.yaml`.

For Step 4d verification: the Personal-View-Customizations variant
must exist. Currently `helpers.playwright.projects.saveCopy({mode='pvc'})`
is NOT registered, so the PVC variant must be created MANUALLY
(via UI Save Copy → mode=Personal View Customizations) before this
manual check, OR via UI as part of this manual session.

## Manual verification steps

### Step 1 — Browse > Dashboards thumbnails preview render quality

1. Browse > Dashboards.
2. Locate the test-project tiles (`AutoTest-CopyClone-<timestamp>`,
   its `-link` companion, and any other variants from chain mode).
3. Verify each tile renders cleanly:
   - **Thumbnail** (dashboard preview OR default-generated thumbnail):
     visible image; NOT a blank box, NOT a broken-image icon, NOT
     visibly mis-sized vs sibling tiles.
   - **Project name** label below thumbnail: visible, no truncation
     artifacts at standard panel widths.
   - **Basic metadata** (creation date / owner / etc.): rendered
     consistently across tiles; aligned with sibling-tile layout.
4. Subjective layout check:
   - All tiles same dimensions, no visual jitter on hover/scroll.
   - No tile overflow / cropping issues.

**Why UI-only:** thumbnail correctness, layout consistency across
tiles, and absence of visual glitches require visual judgment.
Partial DOM probe possible (`img.complete && img.naturalWidth > 0`)
but full render-quality assessment is human.

### Step 4d — View-state customizations preservation on reopen

(Run AFTER the Personal-View-Customizations variant has been created
either manually OR via a future PVC helper.)

1. Open the source project (`AutoTest-CopyClone-<timestamp>`) from
   Dashboards.
2. Apply customizations to its TableView:
   - Add a column filter (e.g., `AGE > 50`).
   - Sort by a column (e.g., `HEIGHT desc`).
   - Hide a column (e.g., toggle off `DIS_POP`).
   - (Optional) Reposition a viewer or change a viewer property.
3. File → Save Project → choose **Personal View Customizations** mode.
   Save as `AutoTest-CopyClone-<timestamp>-pvc`. OK.
4. closeAll.
5. Reopen `AutoTest-CopyClone-<timestamp>-pvc` from Dashboards.
6. Verify view-state preservation:
   - Filter (`AGE > 50`) is active and applied (rows filtered
     appropriately).
   - Sort order (`HEIGHT desc`) is applied (visible row order).
   - Hidden column (`DIS_POP`) is hidden.
   - Viewer positions / properties match the pre-save state.
7. Subjective layout check: panel positions, viewer dimensions, and
   overall workspace layout match the pre-save state. No drift.

**Why UI-only:** filter/sort/column-hide DOM-state IS partially
probable via JS API (`grok.shell.tv.dataFrame.filter.trueCount`,
`tv.grid.columns[i].visible`, sort indicators), but layout positioning
+ viewer dimensions + workspace overall preservation is visual
judgment without an explicit per-element automated criterion.

## Sign-off

PASS — Step 1 thumbnail render + Step 4d customization preservation
both clean; no visual regressions vs prior dev release.
FAIL — list affected step(s) + screenshot; log in
`projects-copy-clone-run.md` under a new dated entry.

## Cross-cutting notes

- Step 1 thumbnail check is best done OPPORTUNISTICALLY when running
  any other manual scenario in this section (deleting / share-project /
  etc.) — covers all visible projects in one pass rather than per-spec.
- Step 4d verification is gated on PVC mode being exercisable. When
  `helpers.playwright.projects.saveCopy({mode})` is registered, PVC
  variant creation can be automated; the visual preservation check
  remains UI-only.

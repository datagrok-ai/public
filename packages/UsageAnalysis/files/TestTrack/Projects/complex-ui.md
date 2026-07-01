---
feature: projects
companion_to: complex.md
companion_spec: (none — complex.md has no single -spec.ts; sibling specs are complex-derived-tables-spec.ts / complex-rename-spec.ts / complex-share-second-user-spec.ts)
ui_only: true
moved_steps_from_canonical: ["Step 4 (drag-drop in Dashboards)", "Step 10 (Move entities to file share / Space)"]
generated_at: 2026-05-03
generated_by: b2-2026-05-03-drag-drop-ui-only-reclassification (Session B Round 2 + Olena verification of Move-to UI option absence)
---

# Complex — UI-only manual companion

This file captures the steps from `complex.md` that cannot currently
be automated through Playwright + JS API mechanisms. Step numbering
preserved from the canonical `complex.md` for cross-reference.

These manual steps are NOT exercised by any of the sibling
`complex-*-spec.ts` files. Human QA must run them as a companion when
validating the complex-mega scenario end-to-end.

## Pre-requisite

The fixture state established by `complex.md` Steps 1–3 (initial
project saved with Data Sync ON, plus 4 additional tables opened from
Spaces / Files / Query / DB) must already be in place when these
manual steps run. Either: (a) drive Steps 1–3 manually first, OR (b)
run `complex-derived-tables-spec.ts` to set up an analogous fixture
state.

## Manual verification steps

### Step 4 — Add newly opened tables to the opened project (drag-and-drop in the Dashboards)

With the project from Step 2 open in the Dashboards view, drag each of
the 4 newly-opened tables (from Step 3) onto the project's Dashboards
tile / project tree node; on drop, the platform should offer
**Link / Clone / Move / Copy** options — pick **Link** (or accept
default) to add as project relations.

**Verify on Context Panel:** the project's relations / tables list
now includes the 4 added tables.

**Why UI-only:** drag-drop event registration is not automatable
through current Playwright + JS API mechanisms. Three mechanisms were
attempted in Session B + Round 2 (synthetic `DragEvent`, raw
`mousedown`/`mousemove`/`mouseup` sequence, CDP-level
`Input.dispatchMouseEvent` via `chrome-devtools__drag` tool) — none
triggered Datagrok's internal drop handlers. Hypothesis: source-side
payload registration via `grok.events` is not exposed to Playwright
without internal-API hooks. See decision-log
`b2-2026-05-03-drag-drop-ui-only-reclassification` and Session B
Round 2 batch report for full forensics.

### Step 10 — Move the following entities to any file share, then to any Space

Each entity is moved twice: first to a file share, then to a Space.

- **Script** — drag-and-drop in Browse tree, OR via right-click
  → **Move to** → select target file share, then move-again to a
  Space.
- **Query** — same flow.
- **Project** — same flow.

**Why UI-only:** both documented paths are unavailable for automation
in the current state:

- **Drag-drop in Browse tree** — same blocker as Step 4 (drag-drop
  registration mechanism not accessible from Playwright; see Step 4
  rationale above).
- **Right-click → Move to** — this menu option does **not exist** in
  the current Datagrok UI (verified by Olena 2026-05-03 against
  dev.datagrok.ai). The original test scenario references a UI
  affordance that has been removed or never existed in the
  right-click context menu.

Therefore Step 10 cannot be automated through any documented path
until either (a) the UI gains a `Move to` menu option in the
right-click context menu, or (b) the drag-drop registration mechanism
is understood and made addressable from Playwright.

The Wave 1b complex-move sub-spec was already DROPPED rather than
inline fragile selector code (see decision-log
`wave-1b-complex-split-b70-followup`). This entry formalizes that
drop and routes the manual coverage here.

## Cleanup responsibility

After the manual run, restore entity locations to baseline so that
downstream chain steps (Step 11 verify, Step 12 share, Step 13
re-auth) operate on the expected entity namespaces:

- Move Script back to its original namespace (`Samples:Cars`).
- Move Query back to its original namespace
  (`Samples:PostgresCustomers`).
- Move Project back to its original Dashboards location.

Log the cleanup actions in `complex-run.md` under the dated entry for
this manual run.

## Sign-off

PASS — both manual steps completed; Step 4 drop-action chooser
appeared and Link selection succeeded; Step 10 entity-move dialogs
operated correctly; baseline restored; no visual regressions.

FAIL — list affected step(s) + screenshot of failure (drop chooser
absent / Move-to menu option absent / dialog error / etc.); log in
`complex-run.md` under a new dated entry. If Step 10's right-click
**Move to** menu option appears in the UI during a future manual run,
note it as a finding so the right-click path can be reclassified as
automatable.

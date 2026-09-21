---
feature: projects
realizes: [views.projects]
target_layer: manual-only
coverage_type: regression
companion_to: complex.md
manual_only_reason: |
  Both steps hinge on drag-and-drop (Dashboards tiles and the Browse tree),
  which cannot be driven by automation on the current build — run them by
  hand.
---

# Complex — UI-only manual companion

This file captures the steps from `complex.md` that cannot currently
be automated. Step numbering preserved from the canonical `complex.md`
for cross-reference.


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

### Step 10 — Move the following entities to any file share, then to any Space

Each entity is moved twice: first to a file share, then to a Space.

- **Script** — to a file share: drag-and-drop in the Browse tree; to a
  Space: right-click → **Move to Space...** and pick the target Space.
- **Query** — same flow.
- **Project** — same flow.

## Cleanup responsibility

After the manual run, restore entity locations to baseline so that
the downstream steps of `complex.md` (Step 11 verify, Step 12 share,
Step 13 re-auth) operate on the expected entity namespaces:

- Move Script back to its original namespace (`Samples:Cars`).
- Move Query back to its original namespace
  (`Samples:PostgresCustomers`).
- Move Project back to its original Dashboards location.

---
{
  "order": 1,
  "datasets": ["System:AppData/Chem/tests/spgi-100.csv", "System:DemoFiles/demog.csv"]
}

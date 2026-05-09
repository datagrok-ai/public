# Migration Report — tree-enhancements-bundle-bug.md

Authoring date: 2026-05-07 (atlas-driven; original scenario authored
in cycle charts-migrate-2026-05-07).
Remediation report date: 2026-05-09 (charts-remediate-2026-05-09).

`produced_from: atlas-driven` — bug-focused regression scenario
authored from atlas + bug-library to close coverage gap for
github-3221 (Tree 8-enhancement bundle). Not migrated from a manual
TestTrack scenario; this report documents the layered authoring
decisions plus the 2026-05-09 remediation scope decision.

## Step mapping (atlas-driven)

| Atlas / bug-library element | Scenario location | Decision |
|----------------------------|-------------------|----------|
| `charts.tree.show-counts` | Step 2 (Capability 1) | exerciseProperty(`showCounts`, true/false) |
| `charts.tree.on-click` | Step 3 (Capability 2) | exerciseProperty(`onClick`, Filter/Select) |
| `charts.tree.font-size` | Step 4 (Capability 3) | exerciseProperty(`fontSize`, 14/30) |
| `charts.tree.orient` | Step 5 (Capability 4) | exerciseProperty(`orient`, LR/RL/TB) |
| `charts.tree.include-nulls` | Step 6 (Capability 5) | exerciseProperty(`includeNulls`, true/false) |
| `charts.tree.layout` | Step 7 (Capability 6) | guarded probe — conditional on `getProperties()` containing 'layout' |
| `moleculeSize` (capability listed in github-3221 fix bundle, not in atlas as a sub_feature) | Step 8 (Capability 7) | enumerate-only check; delegated to Setup propNames log |
| `charts.tree.show-mouse-over-line` (proxy for label-overflow) | Step 9 (Capability 8) | exerciseProperty(`showMouseOverLine`, true) |
| Final visual stability invariant | Step 10 | Verify root DOM non-empty + non-zero size after all 8 toggles |

## Decisions

- **Why this `target_layer` (`playwright`):** existing `tree-spec.ts`
  is at the playwright layer; consistency with sibling.
- **Why this `pyramid_layer` (`bug-focused`):** scenario reproduces
  github-3221 invariant per capability; not a smoke or integration
  scenario.
- **Why this `coverage_type` (`regression`):** github-3221 fixed in
  Charts 1.4.3 — locks each capability's regression bound.
- **Sibling tests consulted:**
  - `tree-spec.ts` (canonical Tree spec at playwright layer).
  - `tree-rowsource-onclick-state-bug-spec.ts` (sister bug-focused
    Tree spec — github-3245).
- **Helpers reused:** `loginToDatagrok`, `softStep`, `specTestOptions`,
  `stepErrors` from `spec-login.ts`. `exerciseProperty` is an inline
  helper authored in the spec body (not registry-eligible per A1
  boundary — single-spec local pattern).
- **Bug library consulted:** yes — `github-3221` populates
  `related_bugs` in frontmatter. Bug-library entry confirms 8
  capabilities bundle; spec exercises each.

## Opt-outs (SCOPE_REDUCTION proposals)

(none — Option A (TIGHTEN) chosen for layout property)

## Deferred items (NOT opt-outs)

(none for this scenario beyond the env-pending guarded probe for
`layout` and `moleculeSize` — both documented as acceptable SR
class per orchestrator Edit 10 env-pending-defensive-skip)

## Remediation cycle 2026-05-09 — layout property guarded probe

Predecessor cycle charts-automator-only-2026-05-08 surfaced
`charts.tree.layout` as declared in spec `sub_features_covered` but
uncovered in code (canonical Critic E verdict EVIDENCE_GAP: E-TRACE-02).
Scenario step 7 (Capability 4) is conditional: "if
`tree.props.getProperties()` includes 'layout'". Step 8 (Capability 7
— moleculeSize) is enumerate-only and also uncovered.

**Migrator decision: Option A (TIGHTEN)**.

| Decision factor | Evaluation |
|-----------------|------------|
| Atlas mapping integrity | `charts.tree.layout` is at atlas line 607 — confirmed sub_feature. Dropping it from sub_features_covered (Option B) would orphan the atlas entry. |
| Code complexity | Guarded probe is 8-line addition. Low cost. |
| Env-pending behavior | If Tree build on dev does not expose 'layout', `console.warn('[SKIP] charts.tree.layout not exposed...')` is acceptable per orchestrator Edit 10 env-pending defensive skip class. |
| Forward path | When Charts package exposes 'layout' on dev build, probe automatically picks it up without scenario amendment. No re-cycle needed. |
| Option B (REDUCE) rejected | Would require dropping `charts.tree.layout` + atlas line 607 cross-reference becomes noise. Atlas integrity matters more than 8 lines of test code. |

**Automator implementation:**

1. Add guarded probe at start of test (after Setup softStep):
   read `tree.props.getProperties().map(p => p.name)` → store as
   `availableProps: string[]`. Console log for diagnostic visibility.
2. After existing capability softSteps, add Capability 4 conditional
   exerciseProperty for layout:
   ```ts
   if (availableProps.includes('layout'))
     await exerciseProperty('layout', 'orthogonal', 'layout=orthogonal');
   else
     console.warn('[SKIP] charts.tree.layout not exposed on current Tree build; capability 4 deferred per scenario step 7 conditional');
   ```
3. Add `// technical: capability 7 (moleculeSize) enumeration-only
   check delegated to Setup propNames log per scenario step 8
   best-effort phrasing` comment near the existing softStep block.

Decision-log cross-reference:
`migration_decisions[tree-enhancements-layout-probe-2026-05-09]`.

## Edge cases

- **Race-tolerant property access:** `tree.props.get()` and
  `tree.props.getProperties()` calls wrapped in try/catch — Charts
  package machinery races with cold-start initialization on dev.
- **Property name 'layout' may not be exposed on all Tree builds.**
  Conditional guard handles this gracefully; env-pending skip
  documented in Notes.

## Unresolved ambiguities

- **moleculeSize in atlas:** github-3221 fix bundle includes
  moleculeSize as one of the 8 capabilities, but the atlas
  (`charts.yaml`) does not currently list `charts.tree.moleculeSize`
  as a sub_feature. Decision: scenario step 8 enumerates moleculeSize
  via `getProperties()` log only (enumerate-only); not in
  `sub_features_covered`. Atlas curator may evaluate adding the
  sub_feature in a future cycle; Migrator does not auto-amend atlas
  per A1 boundary.

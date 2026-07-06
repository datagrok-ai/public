---
feature: biostructureviewer
sub_features_covered:
  - biostructure.api.viewPdbById
  - biostructure.api.viewPdbByData
target_layer: apitest
coverage_type: regression
produced_from: atlas-driven
related_bugs: []
realized_as:
  - js-api-extension-api.ts
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T14:00:00Z
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T14:05:00Z
    failure_keys: []
    review_round: 1
  b:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-automate-01
    timestamp: 2026-06-04T18:42:00Z
    spec_runs:
      - spec: js-api-extension-api.ts
        result: passed
        attempts: 3
        duration_seconds: 18
        failure_keys: []
  e:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-06-04-biostructureviewer-automate-01
    timestamp: 2026-06-04T15:30:00Z
    failure_keys: []
---

# BiostructureViewer — JS API extension (`viewPdbById` / `viewPdbByData`)

Coverage-extension scenario authored under cycle
`2026-06-04-biostructureviewer-migrate-01` Gate F SR routing. Extends the
section beyond the Mol*-centric smoke, the five bug-focused regression guards,
and the four prior coverage extensions (NGL viewer, Property surface +
edge, Context-panel widgets, Mol* viewport overlay buttons) to the two
remaining package-level JS API entry points for opening a structure
programmatically:

- `BiostructureViewer:viewPdbById(pdbId)` — opens a PDB by ID via the
  `byId` helper, wrapping the call in a task-bar progress indicator
  (atlas `biostructure.api.viewPdbById`, source
  `public/packages/BiostructureViewer/src/package.ts#L105`).
- `BiostructureViewer:viewPdbByData(pdbData, name)` — opens a PDB from a
  raw PDB string plus a structure name via the `byData` helper (atlas
  `biostructure.api.viewPdbByData`, source
  `public/packages/BiostructureViewer/src/package.ts#L115`).

The third sibling `viewBiostructure(content, format?, name?)` is already
covered by `property-surface-extension.md` (atlas
`biostructure.api.viewBiostructure`), so this scenario completes the JS API
trio without redundancy.

Net-new contribution against the live covered-union of 49: both ids
`biostructure.api.viewPdbById` and `biostructure.api.viewPdbByData` are
absent from the prior union → **net_new = +2**. With this scenario the
covered-union climbs 49 → 51 = **70.83% / 72**, strictly above the
F-STRUCT-COVERAGE-01 70% threshold — the breadth gap clears on this
round.

`target_layer` rationale (STEP D non-UI tail rule): both functions are
registered via `@func` / `grok.functions.call('BiostructureViewer:…', …)`
in `public/packages/BiostructureViewer/src/package.g.ts` (declarative
`name:` annotations on the package functions). They are grok-callable JS
API entry points whose behaviour is observable through return values and
the Datagrok shell state (`grok.shell.v`, `grok.shell.tv` after the call
docks the viewer). No UI-only assertion is required to validate the API
contract; an `apitest` exercise is the canonical layer for this surface.

`coverage_type` rationale (STEP E): standard JS API exercise of a public
package function trio; not a critical golden path (smoke) and not a
boundary/negative case (edge). Atlas declares no `edge_cases[]` entry
mapping onto these two ids, so the STEP E heuristic applies cleanly →
`regression`.

## Setup

1. Ensure the BiostructureViewer package is loaded (loaded automatically
   on first `grok.functions.call('BiostructureViewer:*', …)` per the
   package-init contract).
2. Close any open table views before each scenario so the docked viewer
   created by the API call has a clean slate to attach to. Use
   `grok.shell.closeAll()` between scenarios.
3. No DataFrame setup required — both functions open a structure-only
   view; they do not require a host table.

## Scenarios

### Scenario 1: `viewPdbById` opens a PDB by RCSB ID and renders

Steps:
1. Call
   `await grok.functions.call('BiostructureViewer:viewPdbById', {pdbId: '1QBS'})`.
   The call returns once the `byId` helper has fetched the structure
   from RCSB and the Mol* viewer has been docked.
2. Resolve the docked viewer: locate the most recently added Biostructure
   viewer by walking `grok.shell.tv.viewers` (or `grok.shell.v` when the
   API opens the structure as a standalone view) and find the entry whose
   `type === 'Biostructure'`.
3. Await render via the viewer's `awaitRendered(timeoutMs)` method
   (atlas `biostructure.viewer` interactions list), e.g.
   `await viewer.awaitRendered(30000)`.

Expected:
- The `grok.functions.call` resolves without throwing.
- A Biostructure (Mol*) viewer is present in shell state after the call
  (the canonical observable side-effect of the function — see the
  task-bar wrapper at `public/packages/BiostructureViewer/src/package.ts#L105`).
- `awaitRendered(30000)` resolves (rather than timing out) — the
  fetched structure parsed and drew successfully.
- No `Parsed object is empty` error is surfaced (the dedicated
  `viewPdbById` entry point wraps the call with a known structure name
  per `byId`, avoiding the bare-`pdb`-prop pitfall documented in atlas
  `edge_cases[5]`).

### Scenario 2: `viewPdbByData` opens a PDB from raw string + name

Steps:
1. Acquire a known PDB text. Two acceptable paths — pick whichever the
   apitest harness already supports without introducing a new fixture:
   - (a) Fetch a small known PDB through `grok.dapi.fetchProxy` (e.g.
     `https://files.rcsb.org/download/1QBS.pdb`) and `.text()` the
     response into a local string.
   - (b) Use a checked-in fixture PDB text from
     `public/packages/BiostructureViewer/tables/` or
     `public/packages/BiostructureViewer/files/` if present; otherwise
     fall back to path (a).
2. Call
   `await grok.functions.call('BiostructureViewer:viewPdbByData', {pdbData: <text>, name: '1QBS'})`.
3. Resolve the docked viewer the same way as Scenario 1 and await render
   with `viewer.awaitRendered(30000)`.

Expected:
- `grok.functions.call` resolves without throwing.
- A Biostructure (Mol*) viewer is present in shell state after the call.
- `awaitRendered(30000)` resolves — the structure parsed from the raw
  string under the supplied `name` argument.
- No `Parsed object is empty` error — the explicit `name` argument is
  the documented safe path that avoids the bare-`pdb` pitfall
  (atlas `edge_cases[5]`, derived from
  `.claude/skills/grok-browser/references/viewers/biostructureviewer.md#L216`).

## Notes

- target_layer rationale: both sub_features are JS API entry points
  registered as `@func`s in `package.g.ts`. STEP D non-UI tail rule
  routes them to `apitest` rather than fabricating a UI scenario.
- atlas entries derived from source:
  `public/packages/BiostructureViewer/src/package.ts#L105` (`viewPdbById`)
  and
  `public/packages/BiostructureViewer/src/package.ts#L115` (`viewPdbByData`).
- Net-new computation: `net_new = {biostructure.api.viewPdbById,
  biostructure.api.viewPdbByData} − live_covered_union − manual_only[]`
  = `{api.viewPdbById, api.viewPdbByData}` (both absent from prior
  union of 49; atlas `manual_only[]` is empty). After authoring the
  covered-union becomes 51 / 72 = 70.83%, clearing the 70%
  F-STRUCT-COVERAGE-01 threshold.
- The third JS API sibling `viewBiostructure(content, format?, name?)`
  (atlas `biostructure.api.viewBiostructure`) is intentionally not
  re-covered here — `property-surface-extension.md` already covers it
  in tandem with the raw-`pdb` pitfall edge case. Avoiding the overlap
  preserves the redundancy posture (no new strict-subset pair
  introduced against any of the ten prior scenarios).
- Helpers: this scenario does not require any helper from
  `helpers-registry.yaml` — the JS API surface is exercised directly
  via `grok.functions.call`. The shell-state walk to resolve the docked
  viewer is a small inline pattern, not a reusable helper.
- The `awaitRendered` settle pattern (atlas `edge_cases[6]`, derived
  from `.claude/skills/grok-browser/references/viewers/biostructureviewer.md#L50`)
  is mandatory before any render assertion — Mol* parse + draw is
  asynchronous and a bare assertion will race the render pipeline.

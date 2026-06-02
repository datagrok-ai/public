---
feature: bio
sub_features_covered:
  - bio.api.get-bio-lib
  - bio.api.get-seq-handler
  - bio.api.get-helm-monomers
  - bio.api.get-seq-helper
  - bio.api.get-monomer-lib-helper
  - bio.lifecycle.init
target_layer: apitest
coverage_type: smoke
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - bio-service-surface-init-api-spec.ts
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T05:10:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T08:30:00Z
    spec_runs:
      - spec: bio-service-surface-init-api-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 61
        failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T09:00:00Z
    review_round: 1
    failure_keys: []
---

# Bio — Service-surface availability (`getSeqHelper` / `getMonomerLibHelper` / `getBioLib` / `getSeqHandler` / `getHelmMonomers`)

Gate F SR-extension scenario (cycle `2026-06-01-bio-migrate-02`)
that realizes atlas critical_path
`bio.cp.bio-service-surface-init` (`priority: p0`,
`derived_from: public/packages/Bio/src/package.ts#L138`). The
critical_path is the highest-priority entry not yet realized on
the live covered union and is the load-bearing entry of the
first-pass set per Critic F's `scope_reduction_proposal`
(closure path, priority order: critical-path adjacency first).

The scenario exercises Bio's full service-surface — the five
`@grok.decorators.func`-registered functions consumer packages
(Helm, Peptides, BiostructureViewer, Dendrogram, HitTriage)
await during their own init:

(a) `Bio:getSeqHelper` → `ISeqHelper` singleton (atlas
    `bio.api.get-seq-helper`, `package.ts#L1678`).
(b) `Bio:getMonomerLibHelper` → `IMonomerLibHelper`
    (`MonomerLibManager`) singleton (atlas
    `bio.api.get-monomer-lib-helper`, `package.ts#L133`).
(c) `Bio:getBioLib` → merged `IMonomerLib` (legacy compat)
    (atlas `bio.api.get-bio-lib`, `package.ts#L175` — class
    method `static getBioLib()` at L176).
(d) `Bio:getSeqHandler` → per-column `ISeqHandler` (notation,
    splitter, region extractor, stats) (atlas
    `bio.api.get-seq-handler`, `package.ts#L180` — class method
    `static getSeqHandler(sequence)` at L181).
(e) `Bio:getHelmMonomers` → list of HELM monomers used in a
    column (atlas `bio.api.get-helm-monomers`,
    `package.ts#L1240` — function name registration at L1241).

All five functions resolve only after `initBio` (atlas
`bio.lifecycle.init`, `package.ts#L138`) finishes — the package
init function registers `_package.completeInit` after
constructing the SeqHelper + MonomerLibManager singletons and
loading monomer sets. A consumer that awaits any of the five
before init completes blocks on the await rather than erroring.
The grok-browser ref doc `bio.md#L552-L566` documents the
init-order invariant explicitly and is the citation anchor for
the test shape used in this scenario.

## Setup

- Authenticate to Datagrok as the test user (apitest harness;
  the five functions are exposed as Datagrok-runtime
  `grok.functions.call('Bio:<name>', ...)` entry points and
  the runtime serializes them after the Bio package init
  resolves).
- Reuse the canonical Bio test fixtures for the per-column
  handler / HELM-monomer surface:
  - `System.AppData/Bio/tests/filter_HELM.csv` — HELM notation;
    used by Scenarios 2 and 3 to exercise `getSeqHandler` and
    `getHelmMonomers` on a `units=helm` column.
  - `System.AppData/Bio/tests/filter_FASTA.csv` — FASTA
    notation; used by Scenario 2 to cross-check
    `getSeqHandler` returns a notation-correct
    `ISeqHandler` instance on a non-HELM column.
- No project save / load: the service-surface contract is
  in-session; no cleanup required beyond closing the opened
  table views at end of run.

## Scenarios

### Scenario 1 — Service-surface five-function resolution after init

Steps:

1. Trigger Bio package init (no-op if already loaded; in an
   apitest run this is implicit because the test harness loads
   Bio before invoking any registered function). The init
   sequence corresponds to atlas `bio.lifecycle.init`
   (`package.ts#L138`).
2. Await
   `await grok.functions.call('Bio:getSeqHelper')` (atlas
   `bio.api.get-seq-helper`). Capture the returned value as
   `seqHelper`.
3. Await
   `await grok.functions.call('Bio:getMonomerLibHelper')`
   (atlas `bio.api.get-monomer-lib-helper`). Capture as
   `monomerLibHelper`.
4. Await
   `await grok.functions.call('Bio:getBioLib')` (atlas
   `bio.api.get-bio-lib`). Capture as `bioLib`.

Expected:

- All four awaits resolve (no timeout, no thrown error, no
  `null` / `undefined` return value).
- `seqHelper` exposes the `ISeqHelper` shape — at minimum a
  callable `getSeqHandler(column)` method.
- `monomerLibHelper` exposes the `IMonomerLibHelper` shape —
  at minimum a callable `getMonomerLib()` / `getBioLib()`
  accessor.
- `bioLib` exposes the `IMonomerLib` shape — at minimum a
  callable `getMonomer(polymerType, symbol)` accessor.
- No error balloon was raised during the four resolution
  steps (`isErrorBallon` returns false).

### Scenario 2 — `getSeqHandler` returns notation-correct per-column handlers (HELM + FASTA)

Steps:

1. Load `System.AppData/Bio/tests/filter_HELM.csv` via
   `grok.data.loadTable(...)` (or equivalent apitest table
   load). Detector classifies the Macromolecule column with
   `units=helm`.
2. Call
   `const handlerHelm = await grok.functions.call(
   'Bio:getSeqHandler', {sequence: df.col('<helm-col>')})`
   (atlas `bio.api.get-seq-handler`, `package.ts#L180`; the
   static method body at L181 delegates to
   `_package.seqHelper.getSeqHandler(sequence)`).
3. Load `System.AppData/Bio/tests/filter_FASTA.csv`. Detector
   classifies the Macromolecule column with `units=fasta`.
4. Call
   `const handlerFasta = await grok.functions.call(
   'Bio:getSeqHandler', {sequence: df.col('<fasta-col>')})`.

Expected:

- Both calls resolve to a non-null per-column `ISeqHandler`
  instance.
- `handlerHelm` reports HELM-shape notation (its
  `notation` / `units` accessor returns `helm`), and
  `handlerHelm.getSplitter()` / `handlerHelm.getRegion()`
  are callable.
- `handlerFasta` reports FASTA-shape notation (notation
  accessor returns `fasta`).
- The two handlers are distinct instances (per-column
  scoping holds).
- No error balloon raised.

### Scenario 3 — `getHelmMonomers` returns the HELM-column monomer list

Steps:

1. Reuse the `filter_HELM.csv` table loaded in Scenario 2 (or
   re-open). Identify the Macromolecule column with
   `units=helm`.
2. Call
   `const monomers = await grok.functions.call(
   'Bio:getHelmMonomers', {sequence: df.col('<helm-col>')})`
   (atlas `bio.api.get-helm-monomers`,
   `package.ts#L1240`; the function registration at L1241
   carries `name: 'Bio: getHelmMonomers'`).
3. Verify the returned monomer list against the HELM column
   contents — every monomer symbol in the column's HELM
   strings should appear in the returned list, and no
   monomer in the returned list should be absent from the
   source HELM strings.

Expected:

- The call resolves to an array (or iterable) of monomer
  identifiers (strings).
- The returned list is non-empty for a non-empty HELM
  column.
- The set of returned monomer symbols matches the set
  observable by parsing the HELM strings in the column
  (round-trip / mapping consistency).
- No error balloon raised.

## Notes

- atlas critical_path realized: `bio.cp.bio-service-surface-init`
  (`priority: p0`, `derived_from: public/packages/Bio/src/package.ts#L138`).
- atlas cross-feature interaction also touched (informational, not
  net-new in `sub_features_covered`):
  `bio.x.bio-service-surface-to-other-packages` (`coverage_type:
  smoke`) and `bio.x.bio-service-surface-init-order`
  (`coverage_type: smoke`) — both list the five service-surface
  ids plus `bio.lifecycle.init` and motivate the smoke-layer
  emphasis on init-order resolution rather than per-function
  behavioural depth.
- target_layer rationale: `apitest` — all five functions are
  `@grok.decorators.func`-registered API entry points consumed
  via `grok.functions.call('Bio:<name>', ...)` with no UI surface
  of their own. The grok-browser ref doc's
  `## Service-surface availability — bio.flow.service-surface`
  section (`bio.md#L552-L566`) documents the JS-API call shape
  but does not describe a UI flow; per STEP D's non-UI-tail
  routing rule, this is an apitest scenario, not a manufactured
  UI filler. (Sibling per-feature behavioural depth — e.g. the
  HELM monomer-fingerprint similarity-score path — already lives
  inside other Bio playwright scenarios; this scenario owns the
  service-surface resolution contract specifically.)
- coverage_type rationale: `smoke` — the scenario verifies the
  "did the five entry points resolve at all" contract that
  consumer packages await during their own init. Behavioural
  depth (monomer-library-version-dependent expansion of
  `getBioLib`, notation-specific edge cases on `getSeqHandler`,
  etc.) is regression-scope and lives in sibling Bio scenarios
  and in the per-package consumer tests (Helm / Peptides /
  BiostructureViewer / Dendrogram / HitTriage). Atlas
  `critical_paths[].priority: p0` maps to scenario
  `coverage_type: smoke` per the skill STEP E heuristic.
- Sub-features covered:
  - `bio.api.get-bio-lib` (`package.ts#L175`, function body
    L176-L177) — Scenario 1 step 4.
  - `bio.api.get-seq-handler` (`package.ts#L180`, function
    body L181-L184) — Scenario 2 steps 2 and 4 (HELM + FASTA
    per-column instances).
  - `bio.api.get-helm-monomers` (`package.ts#L1240`,
    registration at L1241, body `static getHelmMonomers(...)`
    at L1244) — Scenario 3 step 2.
  - `bio.api.get-seq-helper` (`package.ts#L1678`) — Scenario 1
    step 2 (anchor for the SeqHelper-singleton contract;
    already in live union via sibling scenarios such as
    `bio-lifecycle-immunum-wasm.md`).
  - `bio.api.get-monomer-lib-helper` (`package.ts#L133`) —
    Scenario 1 step 3 (anchor for the
    MonomerLibManager-singleton contract; already in live
    union via sibling scenarios such as
    `bio-lifecycle-monomer-library.md`).
  - `bio.lifecycle.init` (`package.ts#L138`) — Scenario 1
    step 1 (anchor for the init-completion precondition;
    already in live union via sibling scenarios such as
    `bio-lifecycle-immunum-wasm.md`).
- Net-new ids vs `live_covered_union` (this cycle's
  authoritative 61-id covered set):
  - `bio.api.get-bio-lib` — net-new.
  - `bio.api.get-seq-handler` — net-new.
  - `bio.api.get-helm-monomers` — net-new.
  - `bio.api.get-seq-helper` — already in union (anchor).
  - `bio.api.get-monomer-lib-helper` — already in union (anchor).
  - `bio.lifecycle.init` — already in union (anchor).
  - net_new = 3 ids; satisfies the SR-loop progress-sensitive
    bound (delta > 0) and the STEP C net-new refusal. Projected
    coverage after this iteration's merge: 62/99 (~62.6%),
    advancing toward the 70% threshold while the remaining
    first-pass scenarios are authored on subsequent iterations.
- Manual-only subset: none of the six covered sub_features
  appear in atlas `manual_only[]` (verified against atlas
  revision 3 `manual_only[]` list — none of `bio.api.*` or
  `bio.lifecycle.init` is flagged manual_only).
- Deferrals: none. All three scenarios are observable inside a
  single apitest run via direct `grok.functions.call(...)`
  invocations; no UI dialog / pixel-precision / docker-eviction
  dependency.
- Bug-context (`related_bugs`): none. No atlas
  `known_issues[]` / `edge_cases[]` entry maps to the
  service-surface family directly; the closest related
  cross-feature interaction is
  `bio.x.bio-service-surface-to-other-packages` (no bugs
  attached) and `bio.x.bio-service-surface-init-order` (no
  bugs attached). `related_bugs: []` per atlas.
- See: `.claude/skills/grok-browser/references/bio.md#L552`
  (`## Service-surface availability — bio.flow.service-surface`
  H2; covers all five service-surface sub_features per the JS
  API call table at L558-L565).
- This scenario covers 6 sub_features (`F-STRUCT-DENSITY-01`
  floor: 2; `F-STRUCT-INTERACTION-01` floor: 3 in a
  multi-sub_feature scenario — satisfied).

---
{
  "order": 23
}

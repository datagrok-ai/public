---
feature: bio
target_layer: apitest
coverage_type: smoke
priority: p0
realizes: [bio.cp.bio-service-surface-init]
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

Covers Bio's service-layer API surface — five functions other
packages (Helm, Peptides, BiostructureViewer, Dendrogram, HitTriage)
call during their own initialization:

- `Bio:getSeqHelper` → the `ISeqHelper` singleton
- `Bio:getMonomerLibHelper` → the `IMonomerLibHelper`
  (`MonomerLibManager`) singleton
- `Bio:getBioLib` → the merged `IMonomerLib` (legacy compatibility)
- `Bio:getSeqHandler` → a per-column `ISeqHandler` (notation,
  splitter, region extractor, stats)
- `Bio:getHelmMonomers` → the list of HELM monomers used in a column

All five only resolve after Bio's package init finishes — consumer
packages that await any of them before init completes just block on
the await rather than erroring. This is the highest-priority
untested surface in the section: if the init-order guarantee breaks,
every consumer package that depends on Bio at startup breaks with it.

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
   Bio before invoking any registered function).
2. Await
   `await grok.functions.call('Bio:getSeqHelper')`. Capture the
   returned value as `seqHelper`.
3. Await
   `await grok.functions.call('Bio:getMonomerLibHelper')`.
   Capture as `monomerLibHelper`.
4. Await
   `await grok.functions.call('Bio:getBioLib')`. Capture as
   `bioLib`.

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
   (the static method body delegates to
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
   (the function registration carries
   `name: 'Bio: getHelmMonomers'`).
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

---
{
  "order": 23
}

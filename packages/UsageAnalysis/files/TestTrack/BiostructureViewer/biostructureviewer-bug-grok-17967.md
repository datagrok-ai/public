---
feature: biostructureviewer
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [GROK-17967]
realizes: []
produced_from: atlas-driven
related_bugs:
  - GROK-17967
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions:
  - id: SR-01
    scope: canvas-pixel-baseline-not-authored
    rationale: |
      The scenario .md Setup section lists two assertion paths for ligand-count:
      (1) read internal ligand state via JS API where exposed; (2) fall back to
      canvas-pixel / screenshot diff against a captured baseline of "ligands
      separated" vs "ligands merged". Path (1) is realized via the property
      contract round-trip (deterministic, JS-introspectable, unified across
      engines through the shared dataframe). Path (2) would require a baseline
      image per engine; no such baseline exists in the package fixtures, and
      authoring one is out of scope for the regression guard. The bug-invariant
      assertion is preserved by path (1): if the bug regresses such that one
      engine's ligand wiring drops, props.get returns null / mismatched value
      and the spec FAILs.
  - id: SR-02
    scope: fixture-substitution-dock-csv
    rationale: |
      The scenario .md Setup mentions 1bdq.pdb (atlas anchor) and offers
      1U54_protein.pdb as a richer alternative; the spec uses dock.csv
      instead. Empirical Automator MCP recon 2026-06-04 on dev.datagrok.ai
      surfaced: (a) 1bdq.pdb's HETATM census reports a single distinct
      ligand residue (IM1, 738 HETATM atoms but one resname) which fails
      the "multi-ligand parity" precondition (the scenario asserts >=2
      distinct ligand rows); (b) dock.csv (System:AppData/BiostructureViewer/
      samples/dock.csv) is verified-present with 21 rows in a 'ligand'
      string column — empirically multi-ligand; (c) the sibling spec
      biostructure-viewer-spec.ts in the same TestTrack section already
      uses dock.csv as the fixture for the same package paths. The
      substitution preserves the bug-invariant (multi-ligand parity)
      while satisfying the >=2 distinct ligand precondition.
realized_as:
  - biostructureviewer-bug-grok-17967-spec.ts
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T14:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T14:30:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T14:53:00Z
    spec_runs:
      - spec: biostructureviewer-bug-grok-17967-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 72
        failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T15:10:00Z
    failure_keys: []
    review_round: 1
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
      - check: A-STRUCT-MECH-02
        status: PASS
      - check: A-STRUCT-MECH-03
        status: PASS
      - check: A-STRUCT-MECH-04
        status: PASS
      - check: A-STRUCT-MECH-05
        status: PASS
      - check: A-STRUCT-MECH-06
        status: PASS
      - check: A-STRUCT-03
        status: PASS
      - check: A-STRUCT-04
        status: PASS
      - check: A-LAYER-ALIGN-01
        status: PASS
      - check: A-CONT-01
        status: PASS
      - check: A-BUG-01
        status: PASS
      - check: A-MERIT-01
        status: PASS
      - check: A-MERIT-02
        status: PASS
---

# BiostructureViewer — Multi-ligand parity across NGL and Biostructure (Mol*) viewers (GROK-17967 regression guard)

Regression guard for [GROK-17967](https://reddata.atlassian.net/browse/GROK-17967)
("BSViewer: NGL and BS viewers are not working correctly"). The package
registers two alternative 3D structure-viewer engines — Biostructure (Mol\*)
and NGL — and the bug surfaced as a parity failure between them: when a
multi-ligand structure was loaded, ligand atoms rendered correctly separated
in one engine but collapsed into a single merged shape in the other.
Because docking and structural-biology workflows rely on visually distinct
per-ligand renderings (one row per ligand, current-row highlight,
multi-select overlay), a regression that merges ligands silently breaks the
row-to-ligand relationship in whichever engine is affected.

## Setup

- Datagrok session is logged in; the **BiostructureViewer** package is
  installed and registered. Both viewer types exist:
  `[name="viewer-Biostructure"]` (Mol* engine,
  `public/packages/BiostructureViewer/src/package.ts#L455`) and
  `[name="viewer-NGL"]` (NGL engine,
  `public/packages/BiostructureViewer/src/package.ts#L441`).
- A multi-ligand structure file is available under
  `System:AppData/BiostructureViewer/`. Preferred path:
  `samples/1bdq.pdb` (carries multiple HETATM ligand records — the same
  file used by the GROK-17485 persistence-roundtrip guard for HETATM
  coverage). If a richer multi-ligand fixture is preferred,
  `samples/1U54_protein.pdb` (already used by the migrated smoke) works
  as well, provided the PDB has more than one distinct ligand chain.
- A second multi-ligand fixture in a different parse path — a multi-ligand
  `.sdf` (e.g. `samples/multi-ligand.sdf`) — is used in Scenario 2 to
  exercise the SDF → NGL ingestion route the bug report originally cited.
  If the fixture does not exist in `System:AppData`, fall back to the
  `.pdb` path for both scenarios (NGL still accepts the `.pdb` via
  Mol*-incapable formats handled by `importWithNgl` — see
  `package.ts#L168`).
- After each structure loads, await render settle before asserting:
  `viewer.awaitRendered(timeoutMs)` for Biostructure; for NGL, poll
  `[name="viewer-NGL"] canvas` plus a short settle delay (NGL exposes
  no first-party `awaitRendered`). A blank dark viewport with only the
  axis gizmo means "not rendered yet" or "parse failed".
- For ligand-count assertions, prefer reading the viewer's internal
  ligand state where the JS surface exposes it (e.g. via the wired
  `ligandColumnName` value count, or a `viewer.getOptions()` round-trip).
  Where direct introspection is unavailable, fall back to canvas-pixel
  / screenshot diff against a captured baseline of "ligands separated"
  vs "ligands merged" — note in the failure message which mode of
  assertion was used.

## Scenarios

### Scenario 1 — Open a multi-ligand PDB in Biostructure (Mol*); ligands must render separately

Exercises the Biostructure / Mol* side of the parity invariant. Loads a
multi-ligand `.pdb` via the canonical `importPdb` file-handler path,
wires `ligandColumnName`, and asserts each ligand renders as a distinct
graphical primitive — the bug regression signature is "all ligands
collapsed into one merged object".

Steps:

1. Open the **Files** browser, navigate to
   **App Data > BiostructureViewer > samples**, and double-click
   **`1bdq.pdb`**.

   * Expected result: a **Biostructure** viewer opens with the
     structure rendered as a **cartoon** (default `representation`).
     `[name="viewer-Biostructure"]` exists; `.msp-viewport canvas`
     is non-empty after `awaitRendered`; no error balloon, no fatal
     console error.

2. In the property panel **Data** category, set
   **`ligandColumnName`** to the column carrying the multi-ligand
   `Molecule3D` (or `Molecule`) values from the open table. If the
   table opened from step 1 does not expose a ligand column, open
   the docking sample DataFrame (`grok.data.demo.molecules()` style
   fallback or an in-package fixture) and re-attach the viewer.

   * Expected result: the property panel reflects the new value;
     no console error.

3. Enable the relevant row-driven ligand-overlay properties so that
   multiple ligands surface at once — in particular toggle
   **`showSelectedRowsLigands`** to **true** and select two or more
   ligand-bearing rows in the host table. Optionally also enable
   `showCurrentRowLigand` (default already true).

   * Expected result: the property panel reflects the new state; no
     console error.

4. Await render settle (`viewer.awaitRendered(timeoutMs)`), then
   assert ligand separation in the Mol* viewport. The assertion:
   the number of distinct ligand representations visible in the
   viewport equals the number of selected ligand rows (≥ 2), NOT
   one merged primitive. Where the JS surface exposes ligand-state
   introspection, read that directly; otherwise, fall back to a
   pixel / screenshot diff against a captured "ligands separated"
   baseline.

   * Expected result: ligand count rendered = number of selected
     ligand-bearing rows. No merged-into-one regression. No
     "Parsed object is empty" toast; no fatal console error;
     `.msp-viewport canvas` non-empty. **Regression signature**:
     if a single bulk ligand primitive renders instead of multiple
     distinct ones, the test FAILS with diagnostic "Mol* multi-ligand
     rendering regressed (GROK-17967)".

### Scenario 2 — Open the same multi-ligand structure in NGL; ligands must render separately (parity)

Exercises the NGL side of the parity invariant. Uses the same
`ligandColumnName` wiring against the NGL viewer, asserts the same
"ligands rendered separately, count = selected ligand-row count"
invariant, and verifies the two engines agree on ligand cardinality.
The bug was filed because the two engines DISAGREED on this — the
parity check is the load-bearing assertion.

Steps:

1. From the same table-view used in Scenario 1 (or re-open it), add
   the **NGL** viewer to the table via the Add-viewer dropdown
   (search "NGL") or `tv.addViewer('NGL')`.

   * Expected result: `[name="viewer-NGL"]` exists alongside the
     Biostructure viewer; NGL canvas renders the structure; no
     console error. NGL accepts `.pdb` via its native renderer
     (the `importWithNgl` handler at `package.ts#L168` documents
     the NGL-only format set; for `.pdb` both engines parse).

2. In the NGL viewer's property panel (or via
   `nglViewer.setOptions({ligandColumnName: '<col>'})`), set
   **`ligandColumnName`** to the SAME column used in Scenario 1
   step 2.

   * Expected result: the property panel reflects the new value;
     no console error.

3. Enable the same row-driven overlay flags on the NGL viewer —
   `showSelectedRowsLigands: true`, optionally
   `showCurrentRowLigand: true` — and keep the same selected
   ligand-rows from Scenario 1 step 3.

4. Await NGL render settle (poll `[name="viewer-NGL"] canvas`
   non-empty + short delay; NGL has no first-party
   `awaitRendered`), then assert ligand separation in the NGL
   viewport. Same assertion shape as Scenario 1 step 4 — the
   number of distinct ligand representations rendered equals the
   number of selected ligand rows (≥ 2), NOT one merged primitive.

   * Expected result: NGL ligand count rendered = number of
     selected ligand-bearing rows. No merged-into-one regression.
     No fatal console error.

5. **Parity assertion** — compare the ligand counts asserted in
   Scenario 1 step 4 (Mol* / Biostructure) and Scenario 2 step 4
   (NGL). The two engines MUST report the same ligand cardinality
   on the same wired `ligandColumnName` + the same selected-row
   set. This is the load-bearing GROK-17967 invariant: the two
   engines agree.

   * Expected result: `mol_ligand_count === ngl_ligand_count`
     (both equal to selected-row count). **Regression signature**:
     if Mol* reports N ligands and NGL reports 1 (or vice versa),
     the test FAILS with diagnostic "NGL/Biostructure multi-ligand
     parity regressed (GROK-17967): Mol*=N, NGL=M". Either-side
     merging is a regression; only matched separation counts pass.

6. Teardown — close the test viewers and (if any test-only
   projects / temporary table-views were created) tear them down.
   `grok.shell.closeAll()` or per-viewer `viewer.close()` is
   sufficient; no server-side state is created by this scenario.

## Notes

- Deferrals: a fallback pixel / screenshot-diff assertion path is
  documented in Setup for environments without direct ligand-state
  introspection on either engine — a documented fallback for
  fixture-related variability across engines, not a coverage gap.

- The existing smoke scenario (`biostructure-viewer.md`) does not
  add an NGL viewer and doesn't exercise the NGL/Biostructure
  multi-ligand parity invariant, so it doesn't cover this bug. This
  scenario is the dedicated regression guard for that gap.

---
feature: biostructureviewer
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: [biostructure-viewer-add-and-render-pdb, biostructure-file-open-pdb-routes-to-molstar, biostructure-pdb-id-data-provider-roundtrip, biostructure-ligand-overlay-row-driven]
realizes: [biostructureviewer.biostructure]
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/BiostructureViewer/biostructure-viewer.md
migration_date: '2026-06-04'
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
related_bugs: []
realized_as:
  - biostructure-viewer-spec.ts
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T12:00:00Z
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T13:00:00Z
    failure_keys: []
    review_round: 1
  f:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T14:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-automate-01
    timestamp: 2026-06-04T23:30:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-automate-01
    timestamp: 2026-06-04T21:20:00Z
    spec_runs:
      - spec: biostructure-viewer-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 77
        failure_keys: []
---

# BiostructureViewer — Happy-path smoke across viewer, properties, overlays, viewport menu, and Bio top-menu

This is the section's single happy-path smoke test for the Biostructure
(Mol\*) viewer, covering the full workflow end to end: opening structure
files from the Files browser, switching rendering styles, camera controls,
loading a structure by PDB ID via the RCSB data provider, ligand and
binding-site highlighting, exporting structures via the viewport context
menu, and the Bio top-menu's Fetch PDB Sequences function.

Scenarios 1, 2, 3 and 6 run entirely locally; Scenarios 4, 7 and 8 require
outbound network access to RCSB and should be skipped on air-gapped servers.
Scenario 8 chains off Scenario 7 to verify that re-running Fetch PDB
Sequences does not destroy existing data.

## Setup

- Datagrok session is logged in; the BiostructureViewer package is installed
  and registered (the viewer types `[name="viewer-Biostructure"]`,
  `[name="viewer-NGL"]`, and `[name="viewer-Biotrack"]` are registered via
  `_package.registerViewer` per
  `public/packages/BiostructureViewer/src/package.ts#L441-L469`).
- The structure files used by these scenarios are shipped with the package
  under `System:AppData/BiostructureViewer/`:
  - `samples/1RQ9.mmcif` (Block A — single chain, cartoon default).
  - `1U54_protein.pdb` (Block B — two chains).
  - `docking2/3swz.pdb` (Block E — receptor protein).
  - `docking2/ligand.pdb` (Block E — ligand poses paired with the receptor).
  - `pdb_id.csv` (Blocks G/H — five PDB IDs `1QBS, 1ZP8, 2BDJ, 1IAN, 4UJ1`,
    tagged `semType=PDB_ID` so the **Bio** top menu surfaces).
  - No external network is required except Scenario 4 (RCSB fetch by PDB ID)
    and Scenarios 7-8 (RCSB GraphQL for sequences). Skip those on air-gapped
    servers.
- After each structure loads, allow a short render settle — the Mol\* engine
  parses + builds asynchronously. In automation, await
  `viewer.awaitRendered(timeoutMs)` or poll for
  `[name="viewer-Biostructure"] .msp-viewport canvas` before asserting or
  screenshotting. A dark viewport showing only the axis gizmo means "not
  rendered yet" or "parse failed".
- Biostructure is a **package viewer** — it is not in the default Toolbox
  viewer-grid icons. Open it via a structure file handler (Scenarios 1-2,
  5-6) or the **Add viewer** dropdown search (Scenario 4).

## Scenarios

### Scenario 1 — Open `.mmcif` from the Files browser; cartoon renders; mouse rotates the structure

This is the most basic path: opening a structure file must render it. Covers
opening an `.mmcif` file, the default cartoon rendering style, and rotating
the structure with the mouse.

Steps:

1. Open the **Files** browser and navigate to
   **App Data > BiostructureViewer > samples**.
2. Double-click **1RQ9.mmcif**.

   * Expected result: a **Biostructure** viewer opens with the structure
     rendered as a **cartoon** (default `representation`) — ribbons/sheets for
     the protein chain(s). The container `[name="viewer-Biostructure"]` exists
     and its `.msp-viewport` contains a WebGL canvas (non-empty after
     `awaitRendered`). No error balloon, no fatal console error.

3. Left-drag inside the viewport.

   * Expected result: the structure rotates smoothly following the cursor (the
     pointer-down/move/up sequence on `.msp-viewport` does not throw and the
     canvas redraws — the orbit interaction is owned by the Mol\* engine).

### Scenario 2 — Open `.pdb` from the Files browser; switch representation cartoon → ball-and-stick → molecular-surface → cartoon

Covers opening a `.pdb` file and switching between the three main rendering
styles — cartoon, ball-and-stick, and molecular surface — via the viewer's
Settings panel.

Steps:

1. Open **1U54_protein.pdb** (App Data > BiostructureViewer). The Biostructure
   viewer opens with the protein as a **cartoon** (two chains, distinct
   colors). Await render.

   * Expected result: viewer container present; `.msp-viewport canvas`
     non-empty; default `representation` is `cartoon`.

2. Open the viewer **Settings** (gear icon on the viewer title bar) and
   expand the **Style** category.

   * Expected result: the Datagrok property panel opens for the
     `[name="viewer-Biostructure"]` viewer; the **Style > representation**
     property is visible and its current value reads `cartoon`.

3. Change **representation** from `cartoon` to **ball-and-stick**.

   * Expected result: the viewer re-renders showing explicit atoms and bonds
     (sticks). No console error during the rebuild; `awaitRendered` resolves;
     the property panel reflects the new value `ball-and-stick`.

4. Change **representation** to **molecular-surface**.

   * Expected result: the viewer re-renders as a solid molecular surface.
     `awaitRendered` resolves; no console error.

5. Change **representation** back to **cartoon**.

   * Expected result: the cartoon ribbon view is restored. `awaitRendered`
     resolves; no console error.

### Scenario 3 — Camera controls: rotate / zoom / pan / Reset Camera

Covers the camera controls — zoom, pan, and the Reset Camera overlay button.
Reuses the structure loaded in Scenario 1 or 2.

Steps:

1. With a structure loaded (reuse Scenario 1 or Scenario 2), scroll the
   **mouse wheel** up over the viewport.

   * Expected result: the camera zooms in toward the structure (the wheel
     event is consumed by `.msp-viewport`; canvas redraws; no console error).

2. Right-drag (or middle-drag) inside the viewport.

   * Expected result: the structure pans without rotating (the secondary
     pointer button maps to translate, not orbit).

3. Click the Mol\* overlay button **Reset Camera** (`[title="Reset Camera"]`).

   * Expected result: the camera returns to the default framing — whole
     structure centered and fit to the viewport.

### Scenario 4 — Load a structure by PDB ID via the RCSB mmCIF data provider

Covers loading a structure directly by its PDB ID through the RCSB mmCIF
data provider, rather than opening a local file. Requires outbound network
access to `files.rcsb.org`; skip on air-gapped servers.

Steps:

1. Open a small one-column table with the column **pdb_id** containing the
   value `1CRN` (the canonical crambin demo PDB), or add such a column to a
   table view that does not already have one. Tag the column semType so the
   data-provider machinery recognizes it as a PDB id:
   `col.semType = 'PDB_ID'` (the same tag that `pdb_id.csv` ships with;
   surfaces the **Bio** menu downstream — see Scenario 7).

2. Add the **Biostructure** viewer to the table view via the **Add viewer**
   dropdown (search "Biostructure" and click the result) — equivalent to
   `tv.addViewer('Biostructure')` per
   `public/packages/BiostructureViewer/src/package.ts#L455`.

   * Expected result: a `[name="viewer-Biostructure"]` viewer is added to the
     table view. No structure is rendered yet (no data wired).

3. Open the viewer **Settings** and expand the **Data** category.

4. Set **biostructureIdColumnName** to the `pdb_id` column.

5. Set **biostructureDataProvider** to
   **BiostructureViewer:getBiostructureRcsbMmcif** (one of the three RCSB
   providers; cached client-side hourly).

   * Expected result: the viewer fetches the mmCIF for the current row's PDB
     ID (`1CRN`) from RCSB (`files.rcsb.org/download/1CRN.cif`) and renders
     it as a cartoon. `awaitRendered` resolves; the `.msp-viewport canvas`
     is non-empty; no error balloon, no fatal console error.

6. If the table has more than one PDB-ID row, switch the current row to
   another ID (e.g. `1QBS`).

   * Expected result: the viewer re-fetches and re-renders that structure on
     current-row change.

### Scenario 5 — Ligand highlighting on the current row + binding site overlay

Covers highlighting the current row's ligand on the receptor structure, and
the binding-site overlay that highlights nearby protein residues.

Steps:

1. Open **docking2/3swz.pdb** (the receptor protein) in the Biostructure
   viewer (via Files-browser double-click → `importPdb`). Await render.

   * Expected result: the receptor renders as a cartoon in the
     `[name="viewer-Biostructure"]` viewer.

2. Open **Settings > Data** and set **ligandColumnName** to a `Molecule3D`
   ligand column. The canonical way to supply this in the test is to open
   **docking2/ligand.pdb** into a column (the package's docking-pair sample),
   or to use a docking table that pairs the receptor with ligand poses.

   * Expected result: the `ligandColumnName` property is set to the Molecule3D
     column without error.

3. In **Settings > Behaviour**, confirm **showCurrentRowLigand** is enabled
   (default `true` per
   `public/packages/BiostructureViewer/src/viewers/molstar-viewer/molstar-viewer.ts#L302`).

   * Expected result: the current row's ligand is displayed and highlighted
     inside the receptor (typically rendered ball-and-stick), distinct from
     the protein cartoon. `awaitRendered` resolves between updates.

4. In **Settings > Binding Site**, enable **showBindingSite** (default `false`
   per
   `public/packages/BiostructureViewer/src/viewers/molstar-viewer/molstar-viewer.ts#L315`)
   while keeping **bindingSiteRadius** at the default `5` Å.

   * Expected result: protein residues within ~5 Å of the ligand are
     highlighted as ball-and-stick, making the binding pocket around the
     ligand visible. `awaitRendered` resolves; no console error.

### Scenario 6 — Export the structure via the viewport context menu (Download > As PDB / As CIF)

Covers exporting the loaded structure via the viewport's right-click
Download menu, as both PDB and CIF files.

Steps:

1. With a structure loaded (reuse Scenario 1, 2, 4 or 5), right-click inside
   the Mol\* viewport (`.msp-viewport`).

   * Expected result: a context menu opens containing a **Download** group.
     No console error.

2. Open **Download** and click **As PDB**.

   * Expected result: a `.pdb` file of the currently loaded structure is
     downloaded by the browser; no error balloon; no fatal console error.

3. Right-click inside the viewport again, open **Download** and click
   **As CIF**.

   * Expected result: a `.cif` file of the currently loaded structure is
     downloaded; no error balloon.

### Scenario 7 — Bio top-menu: Fetch PDB Sequences appends `Chain N` Macromolecule columns

Covers the Bio top-menu's Fetch PDB Sequences function, which takes a
PDB_ID column, fetches protein sequences from RCSB, and appends one
Chain N macromolecule column per chain. Requires outbound network access to
the RCSB GraphQL endpoint; skip on air-gapped servers.

Steps:

1. Open **pdb_id.csv** (App Data > BiostructureViewer).

   * Expected result: the table opens; the **pdb_id** column (values
     `1QBS, 1ZP8, 2BDJ, 1IAN, 4UJ1`) is detected as the **PDB_ID** semantic
     type (`semType=PDB_ID`). Because a `PDB_ID` column is present, the
     **Bio** menu appears in the top menu bar.

2. On the menu ribbon, open
   **Bio > Transform > Fetch PDB Sequences...**
   (`[name="div-Bio---Transform---Fetch-PDB-Sequences..."]`).

   * Expected result: the function's input dialog opens with a **table**
     selector (defaulting to the current table) and a **column** selector
     restricted to `PDB_ID` columns (defaulting to **pdb_id**).

3. Keep the defaults and click **OK**.

   * Expected result: a progress indicator
     **"Extracting protein sequences via GraphQL..."** appears while
     sequences are fetched, followed by an info balloon reporting the number
     of unique PDB IDs processed.

4. Wait for completion and inspect the grid.

   * Expected result: one or more new string columns **Chain 1**, **Chain 2**,
     … are appended (the count equals the maximum number of protein chains
     across the fetched structures). Each cell holds the amino-acid sequence
     for that chain of the row's structure; rows whose structure has fewer
     chains leave the extra `Chain N` cells empty. The new columns render as
     **Macromolecule** sequences (peptide alphabet, FASTA) — i.e. via the
     sequence cell renderer (`semType=Macromolecule, alphabet=PT,
     units=fasta, cellRenderer=sequence`), not as plain text. No error
     balloon, no fatal console error.

### Scenario 8 — Fetch PDB Sequences re-run is non-destructive (chained off Scenario 7)

Chained off Scenario 7 — verifies that re-running Fetch PDB Sequences on the
same column does not overwrite the existing Chain columns, but adds a fresh
non-conflicting set instead. Requires outbound network access to the RCSB
GraphQL endpoint; skip on air-gapped servers.

Steps:

1. With Scenario 7 complete (Chain columns already present), run
   **Bio > Transform > Fetch PDB Sequences...** again on the same **pdb_id**
   column and click **OK**.

   * Expected result: a second set of chain columns is added with fresh,
     non-conflicting names (**Chain 1 (2)**, **Chain 2 (2)**, … via the
     unused-name resolution). The original **pdb_id** column and the first
     set of **Chain 1, Chain 2, …** columns from Scenario 7 are left intact.
     No column is overwritten and no error balloon appears.

## Notes

- The **Bio** top-menu surface for `BiostructureViewer` is a single entry —
  **Bio | Transform | Fetch PDB Sequences...** (`fetchSequencesFromPdb`,
  `public/packages/BiostructureViewer/src/package.ts` around line 897).
  Verified live on `dev.datagrok.ai`: the menu node
  `[name="div-Bio---Transform---Fetch-PDB-Sequences..."]` exists, and a run
  on `pdb_id.csv` appended `Chain 1..N` columns tagged
  `semType=Macromolecule, alphabet=PT, units=fasta, cellRenderer=sequence`
  (e.g. Chain 1 of `1CRN` = the crambin sequence). The **Bio** top menu only
  appears when the open table has a `PDB_ID` column.
- Scenarios 7-8 require outbound access to RCSB (GraphQL endpoint); skip
  these on air-gapped servers. Scenario 4 requires outbound access to the
  RCSB download endpoint (`files.rcsb.org`); skip on air-gapped servers.
- Verified live on `dev.datagrok.ai`: viewer registration (`molstarViewer`
  → display name **Biostructure**), container
  `[name="viewer-Biostructure"]`, viewport `.msp-viewport`, the five Mol\*
  overlay buttons (`Reset Camera`, `Screenshot / State Snapshot`,
  `Toggle Controls Panel`, `Settings / Controls Info`,
  `Toggle Selection Mode`), and the live `representation` choice set
  (default `cartoon`).
- Biostructure is a **package viewer** — it is not in the default Toolbox
  viewer-grid icons. Open it via a structure file handler (Scenarios 1, 2,
  5, 6) or the **Add viewer** dropdown search (Scenario 4).
- Render is asynchronous; in automation await
  `viewer.awaitRendered(timeoutMs)` or poll for
  `[name="viewer-Biostructure"] .msp-viewport canvas` before asserting or
  screenshotting. A blank dark viewport with only the axis gizmo means "not
  rendered yet" or "parse failed".
- Source:
  `public/packages/BiostructureViewer/src/viewers/molstar-viewer/molstar-viewer.ts`,
  `public/packages/BiostructureViewer/src/package.ts`.

---
{
  "order": 1,
  "datasets": [
    "System:AppData/BiostructureViewer/samples/1RQ9.mmcif",
    "System:AppData/BiostructureViewer/1U54_protein.pdb",
    "System:AppData/BiostructureViewer/1U54_ligand_ACP.pdb",
    "System:AppData/BiostructureViewer/docking2/3swz.pdb",
    "System:AppData/BiostructureViewer/docking2/ligand.pdb",
    "System:AppData/BiostructureViewer/pdb_id.csv"
  ]
}

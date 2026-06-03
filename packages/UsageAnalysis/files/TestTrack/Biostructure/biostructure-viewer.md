# Biostructure (Molstar) viewer — happy-path smoke

## Setup

- These scenarios load structure files from `System:AppData/BiostructureViewer/`
  (shipped with the package). No external network is required except Block D
  (RCSB fetch by PDB ID).
- After each structure loads, allow a short render settle (the Mol\* engine
  parses + builds asynchronously); a dark viewport showing only the axis gizmo
  means it has not finished rendering yet.

## Scenarios

### Block A — Open a structure file from the Files browser

1. Open the **Files** browser and navigate to
   **App Data > BiostructureViewer > samples**.
2. Double-click **1RQ9.mmcif**.
* Expected result: a **Biostructure** viewer opens with the structure rendered
  as a **cartoon** (default representation) — ribbons/sheets for the protein
  chain(s). The Mol\* viewport (`.msp-viewport`) contains a WebGL canvas. No
  error balloon, no fatal console error.
3. Left-drag inside the viewport.
* Expected result: the structure rotates smoothly following the cursor.

### Block B — Add the viewer to a table and switch representation

1. Open **1U54_protein.pdb** (App Data > BiostructureViewer). The Biostructure
   viewer opens with the protein as a **cartoon** (two chains, distinct colors).
2. Open the viewer **Settings** (gear icon on the viewer title bar) and expand
   the **Style** category.
3. Change **representation** from `cartoon` to **ball-and-stick**.
* Expected result: the viewer re-renders showing explicit atoms and bonds
  (sticks); no console error during the rebuild.
4. Change **representation** to **molecular-surface**.
* Expected result: the viewer re-renders as a solid molecular surface.
5. Change **representation** back to **cartoon**.
* Expected result: the cartoon ribbon view is restored.

### Block C — Camera controls (rotate / zoom / reset)

1. With a structure loaded (reuse Block A or B), scroll the **mouse wheel** up
   over the viewport.
* Expected result: the camera zooms in toward the structure.
2. Right-drag (or middle-drag) inside the viewport.
* Expected result: the structure pans without rotating.
3. Click the Mol\* overlay button **Reset Camera** (`[title="Reset Camera"]`).
* Expected result: the camera returns to the default framing (whole structure
  centered and fit to the viewport).

### Block D — Load a structure by PDB ID via a data provider

1. Open a small table that has a column of PDB IDs (e.g. a one-column table with
   value `1CRN`), or add such a column. Add the **Biostructure** viewer to it via
   the **Add viewer** dropdown (search "Biostructure").
2. Open the viewer **Settings** → **Data** category.
3. Set **biostructureIdColumnName** to the PDB-ID column.
4. Set **biostructureDataProvider** to
   **BiostructureViewer:getBiostructureRcsbMmcif**.
* Expected result: the viewer fetches the structure for the current row's PDB ID
  from RCSB and renders it as a cartoon. Switching the current row to another ID
  re-fetches and re-renders that structure. No error balloon.

### Block E — Ligand highlighting on the current row (docking)

1. Open **docking2/3swz.pdb** (the receptor protein) in the Biostructure viewer.
2. Open **Settings** → **Data** and set **ligandColumnName** to a `Molecule3D`
   ligand column (e.g. open **docking2/ligand.pdb** into a column, or use a
   docking table that pairs the receptor with ligand poses).
3. In **Settings** → **Behaviour**, confirm **showCurrentRowLigand** is enabled.
* Expected result: the current row's ligand is displayed and highlighted inside
  the receptor (typically ball-and-stick), distinct from the protein cartoon.
4. In **Settings** → **Binding Site**, enable **showBindingSite** (keep
   **bindingSiteRadius** at `5` Å).
* Expected result: protein residues within ~5 Å of the ligand are highlighted
  (ball-and-stick), making the pocket around the ligand visible.

### Block F — Export the structure (context menu)

1. With a structure loaded, right-click inside the Mol\* viewport.
* Expected result: a context menu opens containing a **Download** group.
2. Open **Download** and click **As PDB**.
* Expected result: a `.pdb` file of the current structure is downloaded; no
  error balloon.
3. Right-click again, open **Download** and click **As CIF**.
* Expected result: a `.cif` file of the current structure is downloaded.

### Block G — Bio top menu: Fetch PDB Sequences

The only `BiostructureViewer` function exposed on the **Bio** top menu is
**Bio | Transform | Fetch PDB Sequences...** (`fetchSequencesFromPdb`). It takes a
`PDB_ID` column, fetches protein sequences from RCSB (GraphQL), and appends one
`Chain N` macromolecule column per chain.

1. Open **pdb_id.csv** (App Data > BiostructureViewer).
* Expected result: the table opens; the **pdb_id** column (values like `1QBS`,
  `1ZP8`, `2BDJ`, `1IAN`, `4UJ1`) is detected as the **PDB_ID** semantic type.
  Because a PDB_ID column is present, the **Bio** menu appears in the top menu bar.
2. On the menu ribbon, open **Bio > Transform > Fetch PDB Sequences...**
   (`[name="div-Bio---Transform---Fetch-PDB-Sequences..."]`).
* Expected result: the function's input dialog opens with a **table** selector
  (defaulting to the current table) and a **column** selector restricted to
  `PDB_ID` columns (defaulting to **pdb_id**).
3. Keep the defaults and click **OK**.
* Expected result: a progress indicator **"Extracting protein sequences via
  GraphQL..."** appears while sequences are fetched, followed by an info balloon
  reporting the number of unique PDB IDs processed.
4. Wait for completion and inspect the grid.
* Expected result: one or more new string columns **Chain 1**, **Chain 2**, …
  are appended (the count equals the maximum number of protein chains across the
  fetched structures). Each cell holds the amino-acid sequence for that chain of
  the row's structure; rows whose structure has fewer chains leave the extra
  `Chain N` cells empty. The new columns render as **Macromolecule** sequences
  (peptide alphabet, FASTA) — i.e. via the sequence cell renderer, not as plain
  text. No error balloon, no fatal console error.

### Block H — Fetch PDB Sequences: re-run is non-destructive

1. With Block G complete (Chain columns already present), run **Bio > Transform
   > Fetch PDB Sequences...** again on the same **pdb_id** column and click **OK**.
* Expected result: a second set of chain columns is added with fresh,
  non-conflicting names (**Chain 1 (2)**, **Chain 2 (2)**, … via unused-name
  resolution); the original **pdb_id** column and the first set of Chain columns
  are left intact. No column is overwritten and no error balloon appears.

## Notes

- The Bio top-menu surface for `BiostructureViewer` is a single entry —
  **Bio | Transform | Fetch PDB Sequences...** (`fetchSequencesFromPdb`,
  `package.ts` line ~897). It is verified live on `dev.datagrok.ai`: the menu
  node `[name="div-Bio---Transform---Fetch-PDB-Sequences..."]` exists, and a run
  on `pdb_id.csv` appended `Chain 1..N` columns tagged
  `semType=Macromolecule, alphabet=PT, units=fasta, cellRenderer=sequence`
  (e.g. Chain 1 of `1CRN` = the crambin sequence). The **Bio** top menu only
  appears when the open table has a `PDB_ID` column.
- Blocks G/H require outbound access to RCSB (GraphQL); skip on air-gapped servers.
- Verified live on `dev.datagrok.ai`: viewer registration
  (`molstarViewer` → display name **Biostructure**), container
  `[name="viewer-Biostructure"]`, viewport `.msp-viewport`, the five Mol\*
  overlay buttons (`Reset Camera`, `Screenshot / State Snapshot`,
  `Toggle Controls Panel`, `Settings / Controls Info`, `Toggle Selection Mode`),
  and the live `representation` choice set (default `cartoon`).
- Biostructure is a **package viewer** — it is not in the default Toolbox
  viewer-grid icons. Open it via a structure file handler (Block A/B) or the
  **Add viewer** dropdown search (Block D).
- Render is asynchronous; in automation await `viewer.awaitRendered(timeoutMs)`
  or poll for `[name="viewer-Biostructure"] .msp-viewport canvas` before
  asserting/screenshotting.
- Block D requires outbound access to RCSB; skip it on air-gapped servers.
- Source: `public/packages/BiostructureViewer/src/viewers/molstar-viewer/molstar-viewer.ts`,
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

# Chem | Context Panel — External Database Search Panels (ChEMBL / Chemspace / PubChem / DrugBank)

Manual test of the **external database integrations** that appear under the
**Databases** group of the molecule Context Panel when a molecule cell is
clicked. Covers four providers:

| Provider  | Package      | Backend                                                         | Sub-panels exercised                                              |
|-----------|--------------|----------------------------------------------------------------|------------------------------------------------------------------|
| ChEMBL    | `ChemblAPI`  | EBI ChEMBL web service (`www.ebi.ac.uk/chembl/api/data`)        | Substructure Search API, Similarity Search API                   |
| Chemspace | `Chemspace`  | Chemspace REST API (`api.chem-space.com`, OAuth token)         | Similar, Substructure (under the **Databases \| Chemspace** form) |
| PubChem   | `PubChemApi` | NCBI PUG REST (`pubchem.ncbi.nlm.nih.gov/rest/pug`, async)     | Info, Substructure Search, Similarity Search, Identity Search    |
| DrugBank  | `DrugBank`   | **Bundled** DrugBank structure set (offline, in-package)      | Substructure Search, Similarity Search                           |

All HTTP traffic for the three online providers is routed through the platform
proxy (`/api/connectors/proxy?url=...`) — never a raw browser request. DrugBank
runs entirely client-side against a structure file shipped with the package.

## Preconditions

1. **Server with the four packages published and configured.** ChEMBL and
   PubChem need only outbound network access (no credentials). **Chemspace
   requires a server-side API key/credential** to be configured for the
   `Chemspace` connection — without it the Chemspace panel renders an
   authentication error instead of results. DrugBank needs no configuration.
2. **A user logged in** (any role that can open Demo Files and use the Chem
   package).
3. **Chem package loaded** so molecule cells render (RDKit renderer) and the
   Context Panel **Databases** group is registered.

## Datasets

- `System:DemoFiles/chem/smiles.csv` — primary molecule dataset
  (`canonical_smiles` auto-detects as a Molecule column). Used for the online
  providers, whose chemical spaces are large enough that an arbitrary research
  compound returns hits.
- `System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf` — approved-drug dataset
  used for the DrugBank slice, since DrugBank's curated set returns **No
  matches** for arbitrary research compounds. (Alternatively, sketch a known
  drug such as aspirin `CC(=O)Oc1ccccc1C(=O)O`.)

---

## Scenario 1 — Open the Databases group on a molecule cell

1. Open `System:DemoFiles/chem/smiles.csv`. Wait for the grid to render the
   `canonical_smiles` molecule structures.
2. Click the **first molecule cell** in the `canonical_smiles` column.
   - **Expected:** the right-hand Context Panel switches into *molecule* mode
     and shows the accordion groups **Actions**, **Chemistry**, **Databases**,
     **Biology**, **Structure**.
3. Expand the **Databases** group header.
   - **Expected:** it expands into a nested list of providers:
     **ChEMBL**, **Chemspace**, **PubChem**, **SureChEMBL**, **CDD Vault**,
     **Synthon Search**, **DrugBank**. (This test targets ChEMBL, Chemspace,
     PubChem, DrugBank.)
   - Panes are **lazy-loaded** — no search fires until a provider/sub-panel is
     expanded for the first time.

---

## Scenario 2 — ChEMBL

1. With a molecule cell selected and **Databases** open, expand **ChEMBL**.
   - **Expected:** four sub-panels appear — **Substructure Search API**,
     **Similarity Search API**, **Substructure Search (Internal)**,
     **Similarity Search (Internal)**.
2. Expand **Similarity Search API**.
   - **Expected:** a loader appears, then a grid of up to 20 molecule cards is
     rendered. Each card shows a structure image and a **Score: N.NN** label;
     the first card is the query molecule itself with **Score: 1.00**.
   - **Network (verify if checking):** a request to
     `.../chembl/api/data/similarity/<smiles>/40` goes through
     `/api/connectors/proxy` and returns 200.
3. Expand **Substructure Search API**.
   - **Expected:** a loader, then either a grid of matched structures or
     **No matches** (no error). Request hits
     `.../chembl/api/data/substructure/<smiles>` via the proxy.
4. Hover a result card.
   - **Expected:** tooltip shows the **ChEMBL ID** and "Click to open in ChEMBL
     Database".
5. Click a result card.
   - **Expected:** a new browser tab opens the ChEMBL compound report card
     (`www.ebi.ac.uk/chembl/compound_report_card/<ChEMBL_ID>`).
6. Click the **down-arrow ("Open compounds as table")** icon at the top of the
   results panel.
   - **Expected:** a new table view opens containing the search results
     (named "ChEMBL Similarity Search" / "ChEMBL Substructure Search").

---

## Scenario 3 — Chemspace

> Requires the server-side Chemspace API credential (see Preconditions).

1. Expand **Chemspace** under Databases.
   - **Expected:** a form renders with two inputs — **Ship to country**
     (default *United States*) and **Category** (default *CSCS*) — followed by
     two nested panes: **Similar** and **Substructure**.
2. Expand the **Similar** pane.
   - **Expected:** a loader, then up to ~10 molecule cards each with a
     **Score: N.NN** label.
   - **Network:** an OAuth token request to `api.chem-space.com/auth/token`
     followed by a `POST` search — both via `/api/connectors/proxy`, returning
     200.
3. Expand the **Substructure** pane.
   - **Expected:** loader → grid of substructure matches or **No matches**, no
     error.
4. Change **Category** (e.g. to a different Chemspace catalogue) and/or
   **Ship to country**.
   - **Expected:** the Similar/Substructure result lists refresh for the new
     selection (results re-query; previously expanded panes repopulate).
5. **Negative (if the credential is absent or invalid):** the Chemspace panel
   shows an authentication/error message text instead of crashing the Context
   Panel.

---

## Scenario 4 — PubChem

1. Expand **PubChem** under Databases.
   - **Expected:** nested panes — **Substructure Search**, **Similarity
     Search**, **Identity Search** (and **Info**, which appears only when the
     molecule resolves to a PubChem CID).
2. Expand **Similarity Search**.
   - **Expected:** loader, then molecule cards with **Score** labels (query
     molecule scores 1.00).
   - **Network (async pattern):** the proxy request to
     `.../rest/pug/compound/similarity/smiles/<smiles>/JSON?MaxRecords=20`
     returns **202 Accepted** with a *listkey*; the panel then **polls**
     `.../compound/listkey/<key>/...` until it returns 200. Transient 202/500
     responses during polling are expected and must resolve to results without
     a user-visible error.
3. Expand **Substructure Search**.
   - **Expected:** loader → matched structures or **No matches** (same async
     listkey polling).
4. Expand **Identity Search**.
   - **Expected:** loader → the identical compound if PubChem indexes it,
     otherwise **No matches** — no error either way.
5. (If present) Expand **Info**.
   - **Expected:** PubChem identity/property info for the resolved CID renders.

---

## Scenario 5 — DrugBank (offline / bundled set)

> DrugBank searches a structure set shipped inside the package — **no external
> request is made.** Use an approved-drug molecule so the curated set returns
> hits.

1. Open `System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf` (or sketch aspirin),
   and click a molecule cell for a known drug.
2. Expand **Databases > DrugBank**.
   - **Expected:** two sub-panels — **Substructure Search** and **Similarity
     Search**.
3. Expand **Similarity Search**.
   - **Expected:** loader, then molecule cards with **Score** labels
     (an exact match scores 1.00; similarity threshold is 0.75, up to 20 hits).
4. Expand **Substructure Search**.
   - **Expected:** a grid of DrugBank structures containing the query
     substructure.
5. Hover a result card.
   - **Expected:** tooltip shows the **Common name** and "Click to open in
     DrugBank Online".
6. Click a result card.
   - **Expected:** a new browser tab opens the DrugBank drug page
     (`go.drugbank.com/drugs/<DRUGBANK_ID>`).
7. Click the **"Open compounds as table"** icon.
   - **Expected:** a new table view opens with the DrugBank results.
8. **Negative — research compound:** repeat steps 2–3 on a molecule cell from
   `smiles.csv`.
   - **Expected:** both DrugBank panels render **No matches** gracefully (no
     error) — DrugBank's curated set does not contain arbitrary research
     compounds.

---

## Scenario 6 — Cross-cutting behaviors

1. **Empty / malformed molecule.** Click a blank molecule cell (or a cell whose
   value fails to parse).
   - **Expected:** panels show **SMILES is empty** (empty value) or a malformed
     message ("Molecule string is malformed" for ChEMBL, "Error occurred during
     search. Molecule is possibly malformed" for DrugBank). The Context Panel
     does not crash and other groups stay usable.
2. **Switch molecules.** Click a different molecule cell while the Databases
   panes are expanded.
   - **Expected:** the expanded provider panes re-run their searches for the new
     molecule and replace prior results (no stale results from the previous
     molecule).
3. **Console hygiene.** Throughout the test, the browser console shows no
   uncaught exceptions from the database panels (transient proxy 202/500 for
   PubChem polling are network responses, not JS errors).


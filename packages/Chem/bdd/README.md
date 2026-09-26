# Chem behavioral tests

Gherkin features under `features/`, compiled by `@datagrok-libraries/bdd` (`public/libraries/bdd`)
into the Playwright specs under `generated/` — committed, never edited by hand. They translate the
TestTrack Chem section (`packages/UsageAnalysis/files/TestTrack/Chem/`), one feature per case or
group of cases, each saying in its description what it claims and what it leaves out. 33 features,
all but one `@journey` (the data opened once, the scenarios in order as soft steps).

| Folder           | Features | What it claims |
|------------------|----------|----------------|
| `analyze/`       | activity cliffs (with and without an activity), chemical space, elemental analysis, empty input, MMP, R-groups | each analysis from the top menu on SMILES, V2000 and V3000 molecules: the dialog's defaults, the columns and viewers it adds, the readings the result publishes (`cliffs`, `only cliffs`, the MMP `substitutions` and `pairs`), and an all-empty molecule column handled without an error |
| `calculate/`     | clustering, MPO, properties | BitBIRCH clusters, Cluster MCS and the similarity matrix; MPO Score over a profile; chemical properties, toxicity risks and InChI |
| `filters/`       | substructure card, card state, card interplay, column popup | the substructure card's search types, the sketcher and "Filter as you draw", a card from a cell, after a reset and on a cloned view, with other cards, and a column popup's filter moved to the panel — read through the card's readings on the filter panel (`structure of <column>`, `search type of <column>`) |
| `io/`            | import-export | SDF and MOL2 files opened by their Chem importers and drawn; Save as SDF in V2000 and V3000 and for the filtered rows only, checked in the downloaded file |
| `mpo/`           | profile editor | the MPO Profiles app from the Browse tree: create, edit, save, the pMPO model file, delete |
| `panels/`        | info panels, molecule panels | the Chemistry, Biology, Structure, mixture and highlight panes of a molecule's context panel |
| `projects/`      | chem-state | a saved project brings back the substructure filter and the Scaffold Tree |
| `render/`        | cell actions, molecule rendering | copy as SMILES, molfile V2000/V3000, SMARTS and PNG, export as SVG, and sort by similarity from a molecule cell; molecule, reaction and mixture cells drawn in the grid and in a tooltip |
| `scaffold-tree/` | scaffold tree, functions, colors and limits | building, checking, editing and filtering; coloring, blocked generation with its reason, two tables — through the viewer's `nodes`, `scaffold of node N`, `hits of node N` readings and hit areas |
| `search/`        | similarity, diversity, substructure from the menu | the viewers' properties (metric, fingerprint, limit, size, row source) against the cards they show (`search-results` and `chem-search` readings) |
| `sketcher/`      | cell editor | the sketcher opened from a molecule cell |
| `transform/`     | notation, convert notation once, names to smiles, reactions | the conversions and the columns they add, compared molecule by molecule through RDKit in the page |

One scenario is `@known-failure` (GROK-20956, RDKit's own molblock round trip; it goes with an
RDKit_minimal upgrade); the library's [known-failure audit](../../../libraries/bdd/KNOWN_FAILURES.md)
names the step it stops at and the two defects fixed on 2026-09-23.
Not here: what runs a server-side Python script or container (Curate, Mutate, Butina, Generate
Conformers, Descriptors, Map Identifiers, Synthon Search, the 3D Structure and Gasteiger panes) and
the Identifiers pane's outside lookups — the library's `CLAUDE.md`, "What never becomes a feature".

**Features that draw or type a molecule pin the sketcher** (`the molecule sketcher is
"OpenChemLib"`, the platform's default): the choice is the account's, kept on the server, and an
account that picked Ketcher elsewhere would otherwise change what those features see. The account's
own choice comes back when the feature ends.

What the stand needs: Chem published from the same checkout (the viewers' readings are in its
source), Python scripting (the Scaffold Tree generates its tree with a script), the ChEMBL database
and the Chembl package for Names To Smiles, and EDA for the pMPO training set
(`System:AppData/Eda/drugs-props-train.csv`). The datasets are the package's own files under
`System:AppData/Chem/` and the demo files under `System:DemoFiles/chem/` (`bindings/datasets.ts`).

Run from the package directory against a stand with Chem published. Chem is a package of the
`public/` pnpm workspace, so it resolves the library with nothing to link:

```bash
cd public && grok setup                               # once per checkout: the pnpm workspace
cd libraries/bdd && npm run build                     # the library; dist/ is not committed
cd ../../packages/Chem/bdd
npx grok-bdd run generated --workers=2 --reporter=line
npx grok-bdd run generated/filters                    # one folder
```

The Chem bindings are the package's own screen parts and checks (`bindings/`): RDKit comparisons of
molecule columns (`molecules.ts`), the substructure card's search type and cutoff controls
(`filter-card.ts`), the Scaffold Tree and MMP readiness barriers (`scaffold-tree.ts`), the MPO
profiles on the server (`mpo.ts`) and the dialogs whose OK button has to be caught as it appears
(`dialogs.ts`).

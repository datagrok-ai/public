# Bio behavioral tests

Gherkin features under `features/`, compiled by `@datagrok-libraries/bdd` (`public/libraries/bdd`)
into the Playwright specs under `generated/` — committed, never edited by hand. They replace the
TestTrack Bio specs (`packages/UsageAnalysis/files/TestTrack/Bio/`) and the hand-written specs
under `playwright/`, whose assertions each feature carries or strengthens (the review record is
the `/bdd-translate` skill under `public/.claude/skills/`). The feature count and timings below
are from the last full run recorded here; the 2026-09-22 gap round (every TestTrack case and old
spec re-read against the features) added the FASTA file lifecycle, the project round-trips, the
Monomers view and Match dialog, BILN rendering, the two-column settings of the top menu and the
server side of library deletion — each feature says in its description what it leaves out and
why. The Similarity command used to leave the rows of `filter_HELM` whose length differed from
the reference's empty (GROK-20963, fixed 2026-09-21 in `@datagrok-libraries/bio`;
`calculate/scoring.feature` asserts no blanks). The diversity search's worker chunk failures
of the first round were the distance-matrix service terminating workers it had spawned beyond
the job (fixed in `@datagrok-libraries/ml`).

The folders, by subject (each is a `grok-bdd run generated/<folder>` unit):

| Folder       | Features                                               | What it claims                                                                 |
|--------------|--------------------------------------------------------|--------------------------------------------------------------------------------|
| `analyze/`   | sequence-space, activity-cliffs, msa, composition      | the Bio \| Analyze commands: embeddings and the docked scatter plot with the edited method in its description, cliffs over the Activity column, kalign per cluster, the WebLogo's glyphs as hit areas |
| `transform/` | convert, atomic-level                                  | region, notation (separator, HELM), one Monomer column per position, a V3000 molfile per sequence |
| `calculate/` | scoring                                                | Identity is exactly 1 on the reference row, Similarity peaks there and leaves no cell blank |
| `search/`    | subsequence, similarity, diversity                     | the substructure filter keeps the one row that contains the query (with two sequence columns the dialog's OK does, in `menu/`); the neighbour set changes with the current row; a HELM table gets a HELM subset |
| `annotate/`  | annotate, numbering                                    | the difference column is the two chains joined by `#`, liabilities sit where they say, Manage Annotations lists and drops them, Kabat reaches the aligned column |
| `menu/`      | top-menu                                               | every command under its group, every dialog opens and cancels, the manage views open with their content, the search viewers compute; on a two-column table Subsequence Search filters by the query, Composition binds the chosen column, To Atomic Level runs with Non-linear off |
| `render/`    | renderers, cell-actions, fasta-file                    | the grid reports `helm` / `sequence` / `Monomer` as each column's cell type (BILN included) and paints monomers in their colors, a converted column takes its own notation's renderer (GROK-12164); Copy puts the cell on the clipboard in each of the four notations, the current cell's composition and a monomer's details show on the context panel; a .fasta file opened from the computer is a detected sequence table (GROK-18616), Download > As FASTA... writes it back and the file opens as the same sequences |
| `manage/`    | libraries, collections, monomers                       | the shipped library is loaded; toggling a library checkbox reloads the monomer library and Bio says so (`bio-monomer-lib-loaded`) and leaves the others selected, Add uploads `fixtures/bdd-test-lib.json` through the file chooser and it is still listed when the manager opens again, Delete removes it from every storage; the dialog entry lists what the view lists; a collection is made from the New Collection card, selects on click, is deleted after confirmation; Bio \| Manage \| Monomers is a table of every monomer, Match with Monomer Library offers PEPTIDE/RNA/CHEM |
| `projects/`  | round-trips                                            | Sequence Space's embeddings and scatter plot (GROK-19928), an antibody numbering's aligned column, a HELM table's renderer and library survive a project save and reopen; the numbering run again gives the same result |

From a fresh checkout of `public`, against a local stand on `http://localhost:8888` with this
package published (`npx webpack && grok publish localhost` here — the features need the widget
status the viewers gained on 2026-09-08 and the annotation names):

```bash
cd public/libraries/bdd && npm ci && npm run build   # the library (a path dependency; dist/ is not committed)
npx playwright install chromium                      # its browser, once per machine
cd ../../packages/Bio && npm install                 # the package; npm links the library in and puts grok-bdd in .bin
npx grok-bdd link                                    # ONE Playwright: the library's copy into node_modules (redo after every npm install)
npx grok-bdd run --reporter=list                     # compile --check, then Playwright (one page per worker; --workers N)
npx grok-bdd run generated/analyze                   # one folder
```

The stand needs a platform from `core` at or after 2026-09-08 (`Func.topMenu` in the JS API,
`aria-disabled` on dialog buttons, the grid's `getWidgetStatus`), and the dev key of the
`localhost` entry in `~/.grok/config.yaml`. The datasets are the package's own files under
`System:AppData/Bio/` (`bindings/elements.ts`), published with it. The library features change
the user's library selection and put files on the server, and put both back: the selection is
reset in the Background and at the feature's end, the uploaded library and the created
collection are deleted before and by the scenarios. With two library storages on the stand
(the monomerDomainDB package's and the files) Add asks which one takes the file; the feature
answers "Files" — a stand with the files alone never shows that dialog. PepSeA (a Docker
container), Molecules to HELM (a Python script), the functions called directly (scoring, numbering,
library standardization, the single-sequence converters, the service getters — no UI, so package
tests; the library's `CLAUDE.md`, "What never becomes a feature"), Bio | Folding and the Manage
Monomers view's CRUD are not exercised; projects are
saved through the project API (the ribbon Save dialog and Data Sync are not — see
`projects/round-trips.feature`). `@serial` features (libraries, round-trips) read or toggle the
user's library selection and run one at a time.

Atomic-level conversion selects `HELMCoreLibrary.json` for its standard monomer fixtures and
restores the previous selection afterwards. Custom libraries can redefine symbols such as E
without the R3 attachment required by the branched HELM fixture.

Editing: change a feature, `npx grok-bdd compile`, commit the regenerated spec with it;
`npx grok-bdd list-steps` prints every phrase this package can use, including its own
(`bindings/steps.ts`: the readiness step;
`bindings/annotations.ts`: what a column's annotations and liability hits say;
`bindings/monomer-libs.ts`: the library selection, what the loaded library knows
and where it comes from, the library and collection files on the server;
`bindings/library-files.ts`: the Manage Monomers sketcher's readiness). The cliff count of an activity-cliffs plot is the scatter
plot's own `cliffs` reading (`@datagrok-libraries/ml` publishes it), read with the library's steps.

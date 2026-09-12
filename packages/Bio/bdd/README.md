# Bio behavioral tests

Gherkin features under `features/`, compiled by `@datagrok-libraries/bdd` (`public/libraries/bdd`)
into the Playwright specs under `generated/` — committed, never edited by hand. They replace the
TestTrack Bio specs (`packages/UsageAnalysis/files/TestTrack/Bio/`) and the hand-written specs
under `playwright/`, whose assertions each feature carries or strengthens (the review record is
the `/bdd-translate` skill under `public/.claude/skills/`). Eighteen features, 101 scenarios,
1.9 min in one run (2026-09-08, serial or on four workers alike: four shells booting at once take
30 s each on this stand); one claim is an open product finding, not a test gap: the
Similarity command through its dialog leaves the non-reference rows of `filter_HELM` empty where
the same function called directly scores them (`calculate/scoring.feature`). The diversity
search's worker chunk failures of the first round were the distance-matrix service terminating
workers it had spawned beyond the job (fixed in `@datagrok-libraries/ml`).

The folders, by subject (each is a `grok-bdd run generated/<folder>` unit):

| Folder       | Features                                               | What it claims                                                                 |
|--------------|--------------------------------------------------------|--------------------------------------------------------------------------------|
| `analyze/`   | sequence-space, activity-cliffs, msa, composition      | the Bio \| Analyze commands: embeddings and the docked scatter plot with the edited method in its description, cliffs over the Activity column, kalign per cluster, the WebLogo's glyphs as hit areas |
| `transform/` | convert, atomic-level                                  | region, notation (separator, HELM), one Monomer column per position, V3000 molfiles with no MASS=1 on a heavy atom (GROK-15176) |
| `calculate/` | scoring                                                | Identity is exactly 1 on the reference row, Similarity peaks there, an empty sequence scores to nothing |
| `search/`    | subsequence, similarity, diversity                     | the substructure filter keeps the one row that contains the query; the neighbour set changes with the current row; a HELM table gets a HELM subset |
| `annotate/`  | annotate, numbering                                    | the difference column is the two chains joined by `#`, liabilities sit where they say, Manage Annotations lists and drops them, Kabat reaches the aligned column |
| `menu/`      | top-menu                                               | every command under its group, every dialog opens and cancels, the manage views open, the search viewers compute |
| `service/`   | service-surface                                        | the getters other packages call resolve after init with the methods they use |
| `render/`    | renderers, cell-actions                                | the grid reports `helm` / `sequence` / `Monomer` as each column's cell type and paints monomers in their colors, a converted column takes its own notation's renderer (GROK-12164); Copy as puts the cell on the clipboard in the chosen notation, the current cell's composition and a monomer's details show on the context panel |
| `manage/`    | libraries, collections                                 | toggling a library checkbox reloads the monomer library and Bio says so (`bio-monomer-lib-loaded`), Add uploads `fixtures/bdd-test-lib.json` through the file chooser and Delete removes it; a collection is made from the New Collection card, selects on click, is deleted after confirmation |

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
container), Bio | Folding, the Manage Monomers view's CRUD and project round-trips are not
exercised.

Atomic-level conversion selects `HELMCoreLibrary.json` for its standard monomer fixtures and
restores the previous selection afterwards. Custom libraries can redefine symbols such as E
without the R3 attachment required by the branched HELM fixture.

Editing: change a feature, `npx grok-bdd compile`, commit the regenerated spec with it;
`npx grok-bdd list-steps` prints every phrase this package can use, including its own
(`bindings/steps.ts`: the readiness step, the molfile and alignment readings;
`bindings/annotations.ts`: what a column's annotations and liability hits say;
`bindings/monomer-libs.ts`: the library selection, what the loaded library knows and where it
comes from, the library and collection files on the server).

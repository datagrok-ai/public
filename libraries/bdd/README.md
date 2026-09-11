# @datagrok-libraries/bdd

Behavioral tests for Datagrok packages, written in Gherkin, bound to the u2 and platform vocabulary,
and compiled into Playwright specs that the package owns and commits. The library ships the engine,
the vocabulary and the `grok-bdd` command; a package keeps its features, its own bindings and the
generated specs under `bdd/`.

```gherkin
Feature: MSA workbench
  Background:
    Given user opens the MSA workbench

  Scenario: Running an alignment from the dialog
    When user selects "peptide" in sequence column input in alignment panel
    And user clicks on run msa button in alignment panel
    Then MSA dialog should be visible
    When user clicks on OK button in MSA dialog
    Then status line should contain text "Aligned 5 sequences with kalign"
    And aligned sequences list should have 5 items
```

Nothing in that scenario is a selector. Steps are the verbs, element phrases resolve from the u2
contract (`data-u2`, `data-u2-name`, `data-u2-part`) and from a small registry of names, and the
compiler reports every step, element or dataset it cannot resolve, with the line number.

## Using it in a package

The library is not on npm yet, so a package depends on it by path — what `grok-bdd init` writes
when it runs from a checkout — and npm links the directory into `node_modules`. From a clean clone
of `public`, against a local stand on `http://localhost:8888`:

```bash
cd public/libraries/bdd && npm ci && npm run build   # the library; dist/ is not committed
npx playwright install chromium                      # its browser, once per machine (here, not in the package)
cd ../../packages/<Package> && npm ci                # the package: the library among its dev dependencies
npx grok-bdd link                                    # ONE Playwright: the library's copy into node_modules (redo after every npm ci)
npx grok-bdd run --reporter=list                     # compile --check, then Playwright (4 workers; --workers N, or any Playwright flag)
```

`grok-bdd link` exists because Playwright refuses to be loaded twice in one process and the
library's runtime resolves it from its own directory; the command moves the package's copy to
`node_modules/.bdd-link-backup/` and links the library's in its place (`--undo` puts it back).
Once the library is published, the dependency becomes a version and `link` goes away.

**What the stand needs**: a platform built from `core` at or after 2026-09-10 (the viewer
features rely on signals the core gained for them); a login — global setup mints a token from the
dev key of the `localhost` entry in `~/.grok/config.yaml` (or `DATAGROK_SERVER=<name>`), falls back
to the login form with `DATAGROK_LOGIN`/`DATAGROK_PASSWORD`, and a CI runner passes
`DATAGROK_AUTH_TOKEN`; another stand through `DATAGROK_URL=https://…`. The datasets the features
open are registered in `bindings/platform/datasets.ts`: `demog-1000` is
`System:DemoFiles/demog-1000.csv`, uploaded once with
`grok s files put public/packages/ApiTests/files/datasets/demog-1000.csv "System:DemoFiles/demog-1000.csv" --host localhost`;
`spgi` comes with the published `Chem` package (`grok s packages install Chem`). A sharing feature
shares with the account in `DATAGROK_SHARING_LOGIN`, or, unset, with the `bddsecond` user the
setup creates on the stand (a dev key is needed for that; users cannot be deleted, so it stays).

`grok-bdd init` runs in the package directory and creates what is missing, never overwriting:

```
<package>/bdd/
  package.json        {"type": "module"}   — the specs and bindings here are ES modules
  bdd.config.json     {"tiers": ["viewers"]}   — optional: library tiers beyond the base
  features/**.feature what you write, any hierarchy
  bindings/**.ts      the package's own elements, contexts and steps (optional)
  generated/**.test.ts what `grok-bdd compile` writes — committed, never edited by hand
```

plus `.vscode/settings.json` for the Cucumber extension, the `.gitignore` lines
(`bdd/test-results/`, `bdd/e2e/`, `bdd/.auth.json`), a `test:bdd` script, and `bdd` in the
package tsconfig's `exclude` (webpack type-checks every `.ts` the tsconfig reaches, and `bdd/` is a
Node project with a tsconfig of its own).

```bash
grok-bdd init               # bootstrap bdd/ in the current package (idempotent)
grok-bdd link [--undo]      # wire the package to this checkout of the library
grok-bdd compile            # features/** → generated/**  (--verbose adds how every phrase resolves)
grok-bdd compile --check    # fails when a committed spec is stale — the CI gate
grok-bdd lint               # diagnostics only, notes included
grok-bdd list-steps         # every step this package can use, and where it comes from
grok-bdd run [--headed] [-g "name"] [--reporter=list] [generated/<folder>]
```

Every command runs from the package directory (or from `bdd/`). A feature change needs
`grok-bdd compile` before `grok-bdd run`: the run starts with the drift check and stops on a stale
spec. Results land in `bdd/test-results/` (a trace and a screenshot on failure; `--trace on`
records DOM snapshots too, `--video on` a video); `PLAYWRIGHT_JSON_OUTPUT_NAME=run.json` adds a
JSON report with a duration per step.

## How a feature runs

One browser page per worker: the first scenario the worker runs opens it and boots the shell
(about 4 s), every scenario ends with the shell reset (dialogs and popups closed,
`grok.shell.closeAll()`, the Home view current), and every later feature starts on that reset
shell. A dataset is read from the server once per page and every feature gets a clone, with the
semantic types the first detection found. What a feature leaves on the server it puts back itself
(`atFeatureEnd`). Playwright runs and reports one test per scenario (and per outline row), each
with its own trace; a `Background` runs before every scenario, as Gherkin says.

**`@journey`** on the feature changes that: the feature is one test, the Background runs once, and
the scenarios run in order on the same shell state, each a soft step — a failing scenario is
recorded and the next one still runs, and the test fails at the end listing them. Use it for a
property surface walked section by section, where re-opening the data and the viewer for every
scenario would cost more than the checks. Each scenario then puts back what it changed, and owns
its error and balloon floors. `-g` selects the whole journey.

**`@known-failure`** on a scenario says the product has the defect it describes: its failure does
not fail the test, and its passing does ("the bug is fixed, remove the tag"). Nothing is softened
to stay green.

## Reading a failure

A failed step reports the feature line, the step as written, and the reason in a sentence —
followed, when Playwright gave up waiting for an element, by what the page shows where the phrase
looked. Playwright prints the Gherkin around the line, since the step's location is the feature
file, not the generated spec:

```
StepFailure: features/viewers/box-plot/box-plot.feature:43
  When user right-clicks on the "statsff" area of box plot viewer

Box plot has no "statsff" area right now; it has: view, x axis, y axis, stats, p value, marker

   at ../features/viewers/box-plot/box-plot.feature:43
```

A misspelled menu item gets the visible items, a phrase inside a menu that is not open gets
`context menu: not open`, a misspelled property the nearest captions. A journey lists every failed
scenario in that shape. A stack trace appears only for a programming error in a binding.

## Element phrases

A phrase resolves, in this order, at every level:

| Phrase                                     | Resolution                                                       |
|--------------------------------------------|------------------------------------------------------------------|
| `results`, `browse tab`                    | a registered element (or alias): the whole phrase wins over everything below |
| `second item …`, `last row …`, `3rd input` | an ordinal among the matches                                     |
| `save button in toolbar`                   | composition: `X in|inside|within|on|under Y` — X resolved inside Y (recursively) |
| `label of name input`, `viewers section of toolbox` | `of` names a *part* — of a registered element, or of every element of a kind |
| `sequence column input`                    | a generic **kind** by its longest suffix, qualified by the rest  |
| `"Run MSA" button`, `"First name" input`   | a quoted qualifier: scope words inside it are kept, and a leading "first"/"last" is not read as an ordinal |

**Kinds** cover the whole u2 library (every `data-u2` value it stamps: inputs, buttons, forms,
lists, trees, tables, dialogs, tabs, sections, accordions, cards, wizards, menus, notifications,
tooltips …) plus the Dart shell's conventions (`icon-scatter-plot`, `viewer-Grid`,
`div-section--Viewers`, `input-host-Caption`), and the shell's `viewer`, `view`, `element`. Each
kind knows how its qualifier narrows the candidates — `data-u2-name`, the label or title part,
the text, an aria label, a placeholder, a Dart `name=`. A row goes by its primary text. Every
suffix split is kept: `sequence column input` is tried as `column input` qualified "sequence"
and as `input` qualified "sequence column". `grok-bdd lint` prints how each phrase resolved.

**Look before you phrase.** Open the page in devtools and read the u2 contract — `data-u2`,
`data-u2-name`, `data-u2-part`, the ARIA roles and states — and write the phrase that matches what
is there. A qualifier that starts with an ordinal word needs quotes (`"First name" input`), and a
label equal to a kind name (`Columns input`) resolves to the labelled input first.

**Overrides and contexts.** Register the whole phrase and composition is bypassed:
`workbench.element('save button inside toolbar', {selector: '[data-u2-name="toolbarSave"]'})`.
Element names are global and unique, and the platform base owns the shell's (`toolbox`,
`browse tab`, `context panel`, `grid`, `filter panel`, `context menu` …); a package's own names
live on a *context*, a named region whose vocabulary applies after a step declared with `enters`:

```ts
// bdd/bindings/elements.ts
import {context} from '@datagrok-libraries/bdd';
export const workbench = context('MSA workbench', {selector: '[data-u2-name="msaWorkbench"]'});
workbench.element('results', {selector: '[data-u2-name="results"]'});

// bdd/bindings/steps.ts
import {Given} from '@datagrok-libraries/bdd';
export const openWorkbench = Given('user opens the MSA workbench', async (page) => { … }, {enters: 'MSA workbench'});
```

Inside a context, its names and the generic kinds are looked up in the context root first, then
on the whole page (dialogs and notifications are portaled out of it). A package can also register
a generic kind of its own with `kind()` for a repeated structure that carries no `data-u2`
(U2Demo's `readout`). Selectors inside `labelSelector`, `parts` and `has`-style filters are
evaluated *from the element*: `span:first-child`, not `.x > span:first-child`.

## Steps

`grok-bdd list-steps` prints every phrase the package can use and where it comes from — that
list is the reference; this is the map:

- **Gestures and outcomes on any element** (`bindings/common/steps.ts`): clicks, hovers, typing
  (`types` / `enters` = types and commits), keys, `selects`, checks, expands, drags, `fills in:`;
  `should be/become {state}`, text, value and item counts. States: visible, hidden, present,
  absent, enabled, disabled, checked, unchecked, partially checked, selected, empty, expanded,
  collapsed, focused, invalid, valid — each read from the ARIA state the element uses.
- **The shell** (`bindings/platform/steps.ts`): `user is logged in`, `user opens {dataset}
  dataset` (also `keeping the first N rows [as "name"]`), switching views and table views,
  projects saved and reopened (deleted at feature end), apps, the browse panel, autostarts.
- **The current table through the JS API** (`platform/data.ts`, `columns.ts`): selection and
  filter set and checked row by row, cells, calculated and renamed columns, colour coding,
  other open tables, links between tables, the filter panel's cards through its own API.
- **The top menu and its commands** (`platform/commands.ts`): a path picked by real pointer moves,
  the function call it starts awaited, the columns it added read back.
- **Package functions and their results** (`platform/functions.ts`), **custom platform events**
  (`platform/events.ts`), the clipboard and a file chooser (`common/steps.ts`).
- **The `viewers` tier** (below).

A step definition is an exported `const`; the generated spec imports it by name:

```ts
export const openDataset = Given('user opens {dataset} dataset', async (page: Page, dataset: DatasetEntry) => { … },
  {tier: 'api', description: 'what is checked, when the phrase alone does not say'});
```

Parameter types: `{element}` (any phrase), `{widget}` (a phrase naming a viewer or a widget, the
`grid` or the `filter panel`), `{dataset}`, `{viewer}` (a viewer type), `{key}` (`Enter`,
`Control+A`), `{state}`, plus Cucumber's `{string}`, `{int}`, `{float}`, `{word}`. Runtime helpers
for a package's own steps come from `@datagrok-libraries/bdd/runtime`: `locate`, `gestures.*`,
`viewers.*` (`hitArea`, `hitAreas`, `readValue`, `settle`, `snapshot`, `onViewer`, …),
`expectState`, `expectText`, `atFeatureEnd`.

## Tiers, and the `viewers` tier

`bindings/common` and `bindings/platform` load for every project; vocabulary only some packages
need is a tier under `bindings/tiers/<name>/`, opted into with `{"tiers": ["viewers"]}`.

The `viewers` tier drives viewers the way the platform sees them:

- **Properties by caption** (`user sets "Marker Size Column" property of box plot viewer to
  "WEIGHT"`, `properties of … should be:`), values as text: `true`/`false`, numbers, `#rrggbb`,
  `""` for none, `\n` for a line break.
- **Context menus by path** (`user picks "Misc > Show Inside Values" from the context menu of
  …`, `the open menu should list "Group > Item"`), opened at a named hit area or the `view` area.
- **Hit areas** — the regions a viewer reports in `getWidgetStatus().hitAreas` (`x axis`, `stats`,
  `bar Asian`, `cell 3 of AGE`, `marker`): clicked, double-clicked, hovered, dragged, boxed with a
  key held, typed into, measured, compared.
- **Readings** — `getWidgetStatus().values` (`rows shown`, `bars`, `scene signature`): equal to a
  value, lower/higher/different/the same as before, remembered across a reopen, contained in a
  list reading.
- **Pixels** — `should have repainted [by at least N pixels]`, `less/more ink than before`, an
  area's own ink and repaint, colours in an area (by hue, a shade of anti-aliasing allowed), two
  areas alike or different, the selection highlight (with a margin the selection warrants), the
  value range and the colour scale against before.
- **The legend** (read from its `data-legend-*` attributes), the row tooltip, viewer events,
  layouts saved and loaded, sizes held and restored, and the floors: `no errors should have been
  logged`, `no error or warning balloon should have been shown`.
- **`widgets.ts`** holds the steps first written for one viewer that a second wanted: the viewer's
  own menu, the description's place, empty plot space, range sliders, on-viewer column selectors,
  inner viewers, card readings, lassos, cross-widget drags.

Every property set, menu pick, area gesture, resize and data step snapshots the viewer first, and
no check moves the snapshot, so `should have repainted` and every "than before" compare with the
state before the last change and several checks can follow one change. A settle ends when the
viewer says nothing is pending (`isRenderPending`), not when a cap runs out; a repaint that never
lands is reported as the platform failure it is. **Say what the claim is**: `repainted` is a
change detector, one pixel; a shape gets its own evidence (an area's ink, a colour in an area, a
reading), a chrome toggle takes `by at least N pixels`.

A JS viewer takes part by giving the runtime what a Dart viewer gives it: `getWidgetStatus()`
with its canvas under `parts`, `hitAreas` in CSS px of it and named `values`; a `get
isRenderPending()` true from the render request to the paint; and an `onRendered` observable.

## Generated specs

One `test()` per scenario and per outline row, one `test.step` per Gherkin step located at the
feature line, tags as Playwright tags, `@realizes:<feature>` tags collected into a
`sub_features_covered` header. A `@journey` feature is one `test()` with the Background inline and
every scenario a `run.scenario(...)` soft step. Phrases are emitted as names (`el('…')`), never
selectors, so a selector fix never regenerates anything and `grok-bdd compile --check` fails only
when a feature changed without its spec.

## VS Code

The official Cucumber extension needs, in the package's `.vscode/settings.json`:
`cucumber.features` → `bdd/features/**/*.feature`, `cucumber.glue` → `bdd/bindings/**/*.ts` and
`node_modules/@datagrok-libraries/bdd/bindings/**/*.ts`, and one entry in
`cucumber.parameterTypes` for `state` (`grok-bdd init` writes it from the library's list; the other
custom types are read from the glue). A step shown as undefined while `grok-bdd lint` resolves it
means the settings file is not valid JSON or the glue globs miss the tier directories.

## Developing the library

`npm run build` compiles `src/`, `bindings/` and the Playwright config to `dist/`; `npm run
test:unit` runs the engine tests (nouns, compile, project, init, failure) and the locator test,
which drives the kinds and the platform names over a static page in the library's Chromium (and
skips itself where none is installed); the library is a project itself (`features/platform`):
`npm test` builds, drift-checks it and runs the unit tests — what the libraries CI runs on Node 18,
which is why `@playwright/test` is pinned to the last minor that runs there — and `npm run
test:suite` runs it against a stand. Translating a hand-written spec into a feature, and proving
the feature tests what it claims, is the `/bdd-translate` skill.

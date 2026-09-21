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
cd public && grok setup                              # once per checkout: the pnpm workspace (needs `npm i -g datagrok-tools`)
cd libraries/bdd && npm run build                    # the library; dist/ is not committed (2 s)
npx playwright install chromium                      # its browser, once per machine
cd ../../packages/<Package>/bdd                      # the package's bdd project; the workspace already links the library
npx grok-bdd run --reporter=list                     # compile --check, then Playwright (4 workers; --workers N, or any Playwright flag)
```

Under the pnpm workspace the package and the library resolve to one `@playwright/test`, so
`grok-bdd link` is only for a package installed outside the workspace with npm.

`grok-bdd link` exists because Playwright refuses to be loaded twice in one process and the
library's runtime resolves it from its own directory; the command moves the package's copy to
`node_modules/.bdd-link-backup/` and links the library's in its place (`--undo` puts it back).
Once the library is published, the dependency becomes a version and `link` goes away.

If `grok-bdd run` reports **`Requiring @playwright/test second time`**, run
`npx grok-bdd link` from the package directory and retry. Binding discovery imports both
the library's bindings and the package's bindings, so this can fail during the initial
compile check, before Playwright starts a browser. Two installations of the **same version**
also trigger the error: they must resolve to the same physical copy. Repeat the link after
`npm ci` or another dependency install replaces it; reinstalling the same version alone
does not fix it.

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
records DOM snapshots too, `--video on` a video); `PLAYWRIGHT_JSON_OUTPUT_NAME=<absolute path>.json` adds a
JSON report with a duration per step (a relative name lands in this library's `dist/`, the Playwright
config directory); `node tool/step-times.cjs <bdd project dir> run1.json [run2.json ...]` aggregates
the reports by step phrase and binding — count, total, mean, p90, per run — into `step-times.md`
beside them.

## How a feature runs

One browser page per worker: the first scenario the worker runs opens it and boots the shell
(about 4 s), every scenario ends with the shell reset (dialogs and popups closed,
`grok.shell.closeAll()`, the Home view current), and every later feature starts on that reset
shell. A dataset is read from the server once per page and every feature gets a clone, with the
semantic types the first detection found. What a feature leaves on the server it puts back itself
(`atFeatureEnd`). Playwright runs and reports one test per scenario (and per outline row), each
with its own trace; a `Background` runs before every scenario, as Gherkin says.

Server fixtures can use `{run}` in their names, for example `BDD-Share-Model-{run}`. The suffix is
unique per feature instance (including each worker and repeat) and stays the same across its
scenarios. `{time}` is the feature instance's start in epoch milliseconds, for a name a login must
accept (`[a-z0-9._-]`) and a reader can sort: the users a creation feature makes cannot be deleted,
so they are named by it. String arguments, element phrases, data tables and doc strings resolve both
at runtime;
generated specs stay deterministic. Cleanup registered with `atFeatureEnd` attempts every callback
and fails the run if any callback fails; a run that was killed never gets there, so the server
steps that make or clear a `{run}`- or `{time}`-named fixture, the project save included, also
delete the ones of the same family left by any run older than an hour.

**`@journey`** on the feature changes that: the feature is one test, the Background runs once, and
the scenarios run in order on the same shell state, each a soft step — a failing scenario is
recorded and the next one still runs, and the test fails at the end listing them. Use it for a
property surface walked section by section, where re-opening the data and the viewer for every
scenario would cost more than the checks. Each scenario then puts back what it changed, and owns
its error and balloon floors. `-g` selects the whole journey.

**`@serial`** on the feature runs it one at a time with every other `@serial` feature, while the
rest of the run stays parallel. Use it where features read what other features change at the same
time — a fuzzy gallery search that brings up the fixtures other features create and delete.

**`@known-failure`** on a scenario says the product has the defect it describes: its failure does
not fail the test, and its passing does ("the bug is fixed, remove the tag"). Nothing is softened
to stay green. Outside a journey the tag works on a scenario and on an outline's `Examples` block:
the Background runs plainly, and only the scenario's own steps are the expected failure.

The [known-failure audit](KNOWN_FAILURES.md) records the reproduced defects and the stale tag
removed in September 2026. Inspect the failing step inside each tagged scenario: a green journey
alone does not establish that it failed for the intended reason.

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
| `save button in toolbar`                   | composition: `X in\|inside\|within\|on\|under Y` — X resolved inside Y (recursively) |
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
  collapsed, focused, invalid, valid, ready — each read from the ARIA state the element uses.
  `ready` requires explicit `aria-busy="false"` and no `aria-invalid="true"`; absent readiness
  markup never counts as a completed result.
- **The shell** (`bindings/platform/steps.ts`): `user is logged in`, `user opens {dataset}
  dataset` (also `keeping the first N rows [as "name"]`), switching views and table views,
  projects saved and reopened (deleted at feature end), apps, the browse panel, autostarts; the
  server's spaces, models, groups and roles by name (deleted at feature end), a fixture user (made
  once per stand: users cannot be deleted), its status and who is a (plain or admin) member of a
  group or holds a role; the gallery's render mode and
  its counter against a remembered one (lower, not lower, higher — search, then clear). The membership editor behind Groups..., Roles..., Members
  and Assigned to is `"<name>" membership row` / `membership candidate` with `add button`,
  `remove button` and `checkbox` parts, typed into through `membership search`.
- **The current table through the JS API** (`platform/data.ts`, `columns.ts`): selection and
  filter set and checked row by row, cells, calculated and renamed columns, colour coding,
  other open tables, links between tables, the filter panel's cards through its own API.
- **The top menu and its commands** (`platform/commands.ts`): a path picked by real pointer moves,
  the function call it starts awaited, the columns it added read back.
- **Package functions and their results** (`platform/functions.ts`), **custom platform events**
  and the task bar's progress entries (`platform/events.ts`), the clipboard and a file chooser
  (`common/steps.ts`).
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

The shared gesture helpers interpret `Control` / `Ctrl` as the platform's primary modifier:
Control on Windows/Linux, Command on macOS. Existing `user presses Control+C` and
`… holding Control` steps therefore work on both. This also covers modifier keys held by
drag and legend helpers, and select-all inside the typing and clearing helpers. Internally,
Playwright's `ControlOrMeta` resolves the modifier (requires Playwright 1.45+).

`Delete` / `Del` follows Datagrok's command convention: Backspace on macOS, Delete on
Windows/Linux, so `user presses Shift+Delete` removes selected rows on either platform.
For a physical key, use `ControlLeft` / `ControlRight`, `ForwardDelete` or `Backspace`.
The grid's custom current-cell copy specifically checks physical Control, so that one step
uses `user presses ControlLeft+Shift+C`. `Meta` / `Cmd` and explicit `ControlOrMeta` also work.

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
- **The legend** (read from its `data-legend-*` attributes: its mode — docked, in a corner,
  collapsed to the mini icon, placed nowhere — its slot, its items and their colors, its size
  against before after a splitter drag, whether its items are drawn as structures or as text), the row tooltip, viewer events,
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

Area hovers also settle before the next step. Grid cell tooltip requests participate in the
core viewer's pending-work signal, including the nested correlation grid. Tooltip text checks
consider visible tooltips only; hidden retained text and an absent tooltip satisfy a negative
check. A served core must include the tracked grid tooltip debounce for these absence checks.

A JS viewer takes part by giving the runtime what a Dart viewer gives it: `getWidgetStatus()`
with its canvas under `parts`, `hitAreas` in CSS px of it and named `values`; a `get
isRenderPending()` true from the render request to the paint; and an `onRendered` observable.

A package that customizes an existing widget can contribute its own live areas and readings
through `DG.Widget.addStatusProvider(name, provider)`. For example, Peptides adds the glyphs it
draws in the native grid's headers:

```ts
grid.addStatusProvider('peptides-weblogo', () => ({
  hitAreas: currentGlyphBounds,
  values: {'highlighted rows': highlightedRowCount},
}));
grid.removeStatusProvider('peptides-weblogo');
```

Names identify the contributing component. Registering the same name replaces that provider in
place; later providers override earlier entries. Providers run on each status read and contribute
`parts`, `hitAreas`, and `values`; a JS viewer that overrides `getWidgetStatus()` composes with
`super.getWidgetStatus()`. Geometry uses the widget's coordinate system and must exclude anything
no longer drawn. Detach removes providers; a viewer reattached to a different table needs its
table-specific providers registered again.

## Guides: a scenario as a how-to video

A scenario that answers "how do I …" is also its own demonstration. `grok-bdd guide
features/guides/<name>.feature` compiles it, runs it on one worker in guide mode and renders each
scenario into `guides/<feature slug>/<scenario slug>/`:

- `guide.mp4` — the pointer travels to every element a step acts on, the element is lit (the rest
  of the page dimmed), a small target is zoomed into in place while the click lands, the page
  after the step is revealed, and a caption reads the step as an instruction ("Click on Open local
  file icon in browse toolbar"); a `Then` step shows what it checked with a check mark;
- `step-NN.png` — the lit picture of every step, and `steps.md` — the numbered steps with those
  pictures, ready to paste into a reply;
- with `--gif` also `guide.gif` and `guide-thumb.png`, the docs' own pair.

Guide mode (`BDD_GUIDE=<dir>`, set by the command) records at the step: the page before and after
it (`BDD_GUIDE_SETTLE`, 500 ms by default, lets a dialog or a balloon finish appearing), the last
element the step located, and where the page's own mouse went. A move with a button held is a
drag: the page is pictured along the way (`NN-dragK.png`, at most eight per step), so the video
shows what the drag draws — a selection box, an annotation region — growing under the pointer,
and the step's still shows it complete at the release point. Tests know nothing of it: without
the variable no line of it runs. The viewport is 1080p (1920×1080, so the top menu keeps every group on the bar; `BDD_GUIDE_VIEWPORT=<w>x<h>`
for another) so the video reads without zooming every step, and the shell is the full one (simple
mode off, which the login step and the panel steps read from `shellSimpleMode()`), as a person
has it. Every step is in the video except the login (`guide.silent`) and a step that neither
acted nor changed the page (its before and after pictures are the same file): a table opened
through the API is shown under its caption. A path walked inside a step — the top menu's group,
then each item; a context menu's groups — is a list of stops (`guide.hop`: the page as it was
then, the element's box), and the pointer travels to each with the stop lit; a runtime path that
does not report its stops jumps from the closed menu to the result, so a new menu walk calls
`hop` where `pickTopMenu` and `pickMenuPath` do. Rendering is `tool/guide-render.py` (Pillow + ffmpeg: `py -m
pip install pillow imageio-ffmpeg`, or `FFMPEG=<path>`); `--fps`, `--hold`, `--travel`, `--zoom`
tune the pace when run by hand on a `steps.json` directory.

A feature tagged `@help:<page dir>` (`@help:access/files`) illustrates a help page:
`grok-bdd guide --help-pages` films every such feature and copies each scenario's GIF and thumb
into `<public>/help/<page dir>/img/<scenario slug>.gif` (`BDD_HELP_ROOT` names another tree), so
the walkthroughs on the docs site are regenerated from features rather than recorded by hand.
Guides live under `features/guides/` and run with the rest of the suite: an answer that stops
being true fails a test.

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

Settings are relative to the folder opened in VS Code; nested `.vscode/settings.json` files
are not inherited. With the core repository open, use `public/packages/*/bdd/features/**/*.feature`
for features, and `public/libraries/bdd/bindings/**/*.ts` plus
`public/packages/*/bdd/bindings/**/*.ts` for glue in the root `.vscode/settings.json`.
Include the same `state` parameter type there. With `public/` open, omit the `public/` prefix.

Every settings file that names Cucumber globs also carries `"search.followSymlinks": false`
(`grok-bdd init` writes it). The extension re-scans its globs on every file change through VS
Code's file search, which follows symlinks by default: the pnpm `node_modules` forests and the
junctions `grok-bdd link` makes are circular, so a scan never ends, and they pile up by the dozen
(64 `rg.exe` at three quarters of a 32-core machine, found 2026-09-21). The setting takes effect
on Reload Window; `taskkill /F /IM rg.exe` (or `pkill rg`) clears the ones already running.

## Developing the library

`npm run build` compiles `src/`, `bindings/` and the Playwright config to `dist/`; `npm run
test:unit` runs the engine tests (nouns, compile, project, init, failure) and the locator test,
which drives the kinds and the platform names over a static page in the library's Chromium (and
skips itself where none is installed); the library is a project itself (`features/platform`):
`npm test` builds, drift-checks it and runs the unit tests — what the libraries CI runs on Node 18,
which is why `@playwright/test` is pinned to the last minor that runs there — and `npm run
test:suite` runs it against a stand. Translating a hand-written spec into a feature, and proving
the feature tests what it claims, is the `/bdd-translate` skill.

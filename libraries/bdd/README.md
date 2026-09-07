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

```bash
npm i -g @datagrok-libraries/bdd   # the `grok-bdd` command (until it is on npm: see "From a fresh checkout")
cd <package>
grok-bdd init                      # bootstraps bdd/ — see below — and the manifest entries
npm i                              # installs @datagrok-libraries/bdd and @playwright/test as dev dependencies
grok-bdd run                       # the smoke feature: logged in, Browse visible
```

### From a fresh checkout

The library is not on npm yet, so the packages that have features today (`UsageAnalysis`,
`U2Demo`) depend on it by path — `"@datagrok-libraries/bdd": "file:../../libraries/bdd"`, what
`grok-bdd init` writes when it runs from a checkout of the library — and npm links the directory
into `node_modules` and puts `grok-bdd` in `node_modules/.bin`; `npm ci` needs no registry for
it. Every step, from a clean clone of `public` on any OS — verified 2026-09-07 on Windows:

```bash
# 1. the library: its dependencies, its dist/ (not committed), its browser
cd public/libraries/bdd
npm ci
npm run build
npx playwright install chromium    # the browser of the library's Playwright, once per machine

# 2. the package: its dependencies (the library among them), then ONE Playwright
cd ../../packages/UsageAnalysis    # or U2Demo
npm ci
npx grok-bdd link                  # node_modules/@playwright/test → the library's copy (see below)

# 3. the run
npx grok-bdd compile --check       # the committed specs match the features
npx grok-bdd run --reporter=list   # or: npm run test:bdd
```

`grok-bdd link` exists because Playwright refuses to be loaded twice in one process, and the
library's runtime resolves `@playwright/test` from its own directory (`libraries/bdd/node_modules`)
while the spec resolves it from the package's. The command moves the package's own copy to
`node_modules/.bdd-link-backup/` and links the library's in its place; it is idempotent,
`npx grok-bdd link --undo` puts the copy back, and an `npm ci` in the package undoes it too (run
`link` again after one; without it the run stops at "Requiring @playwright/test second time"). The
browser is installed from the library directory for the same reason: `npx playwright` in the
package would install the build for the package's own Playwright. Once the library is published,
the dependency becomes a version, npm shares the peer `@playwright/test`, and `link` goes away.

**What the stand needs** (the default is `http://localhost:8888`, another stand through
`DATAGROK_URL=https://… npx grok-bdd run`):

- A platform built from `core` at or after `6983855e91` (2026-09-07): the viewer features rely on
  signals the core gained for them — `aria-disabled` on menu items and property rows, the box
  plot's `getWidgetStatus().hitAreas`, single-dash selector names, a menu group opening on the
  first pointer move, and the semantic-type detection changes. On an older platform the
  `viewers` features fail on those steps; the u2 features do not care.
- A login. Global setup mints a token from the dev key of the `localhost` entry (or
  `DATAGROK_SERVER=<name>`) in `~/.grok/config.yaml` — the file `grok config` writes, so a
  `grok publish` setup already has it. Without a key it signs in through the login form with
  `DATAGROK_LOGIN` / `DATAGROK_PASSWORD` (`admin` / `admin` by default), and a CI runner passes
  `DATAGROK_AUTH_TOKEN` directly.
- The datasets the features open (`bindings/platform/datasets.ts`): `demog-1000` is
  `System:DemoFiles/demog-1000.csv` — not part of the demo files, upload it once:
  `grok s files put public/packages/ApiTests/files/datasets/demog-1000.csv "System:DemoFiles/demog-1000.csv" --host localhost`;
  `spgi` is `System:AppData/Chem/tests/spgi-100.csv`, which the published `Chem` package carries
  (`grok s packages install Chem`). Chem also owns the SMILES detector and molecule renderer the
  table-switching scenario exercises.
- For `U2Demo`, the package itself published to the stand (`grok publish localhost` from
  `packages/U2Demo`): its features open the U2 Demo app. `UsageAnalysis`'s viewer features use
  only the platform and the two datasets.

`grok-bdd init` runs in the package directory (it refuses anywhere without a `package.json`) and
creates what is missing, never overwriting what is there: `bdd/package.json`, `bdd/bdd.config.json`,
`bdd/tsconfig.json`, `bdd/bindings/elements.ts` (a context named after the package), `bdd/bindings/steps.ts`
(a `user opens the <Package> app` step to adjust), `bdd/features/smoke.feature`, `.vscode/settings.json` for
the Cucumber extension (merged into an existing one), the `bdd/test-results/`, `bdd/e2e/`, `bdd/.auth.json`
lines in `.gitignore`, and in `package.json` the two dev dependencies and a `test:bdd` script. When the
library is already installed it compiles the smoke feature right away. Running it again reports
everything as existing and changes nothing.

```
<package>/bdd/
  package.json        {"type": "module"}   — the specs and bindings here are ES modules
  bdd.config.json     {"tiers": ["viewers"]}   — optional: library tiers beyond the base
  features/**.feature what you write, any hierarchy
  bindings/**.ts      the package's own elements, contexts and steps (optional)
  generated/**.test.ts what `grok-bdd compile` writes — committed, never edited by hand
```

```bash
grok-bdd init               # bootstrap bdd/ in the current package (idempotent)
grok-bdd link [--undo]      # wire the package to this checkout of the library (until it is on npm)
grok-bdd compile            # features/** → generated/**  (+ errors; --verbose adds how every phrase resolves)
grok-bdd compile --check    # fails when a committed spec is stale — the CI gate
grok-bdd lint               # diagnostics only, notes included
grok-bdd list-steps         # every step this package can use, and where it comes from
grok-bdd run [--headed] [-g "name"] [--reporter=list]   # compile --check, then Playwright
```

Every command runs from the package directory (or from `bdd/` itself). A feature change needs
`grok-bdd compile` before `grok-bdd run`: the run starts with the drift check and stops on a stale spec.

`grok-bdd run` uses the library's Playwright config: `DATAGROK_URL` (default
`http://localhost:8888`), login from `DATAGROK_AUTH_TOKEN` when `grok test` provides it, otherwise
from the developer key of `DATAGROK_SERVER` (default `localhost`) in `~/.grok/config.yaml`, otherwise
the login form with `DATAGROK_LOGIN`/`DATAGROK_PASSWORD`. Results land in `bdd/test-results/`
(traces and screenshots on failure). Add `bdd/test-results/`, `bdd/e2e/` and `bdd/.auth.json` to
the package `.gitignore`. The sample lives in `packages/U2Demo/bdd`: the MSA workbench features at
the root, and under `features/demo/**` one feature per sub-demo of the U2 Demo app — every u2
control the library has a kind for, driven the way its page presents it.

## How a feature runs

One browser page per feature file: the first scenario opens it, every scenario ends with the shell
reset (dialogs and popups closed, `grok.shell.closeAll()`, the Home view current), the last one
closes it. Playwright still runs and reports one test per scenario (and per outline row), so
`-g`, tags, retries, traces and screenshots work as usual. A `Background` runs before every
scenario, as Gherkin says; `user is logged in` only navigates when the page is not in the shell yet.

**`@journey`** on the feature changes that: the feature is one test, the Background runs once, and
the scenarios run in order on the same shell state, each a soft step — a failing scenario is
recorded and the next one still runs, and the test fails at the end listing them. Use it for a
property surface walked section by section, where re-opening the data and the viewer thirteen
times would cost more than the checks (the box plot: 13 scenarios in 14 s, 6 of them the login and
the open). Each scenario then puts back what it changed, so the next starts where the Background
left off. `-g` selects the whole journey; the report shows every scenario and step under it.

## Reading a failure

A failed step reports the feature line, the step as written, and the reason in a sentence —
followed, when Playwright gave up waiting for an element, by what the page shows where the phrase
looked. Playwright prints the Gherkin around the line, since the step's location is the feature
file, not the generated spec:

```
StepFailure: features/viewers/box-plot.feature:43
  When user right-clicks on the "statsff" area of box plot viewer

Box plot has no "statsff" area right now; it has: view, x axis, y axis, stats, p value, marker

   at ../features/viewers/box-plot.feature:43
   42 |       | Show P Value          | true  |
 > 43 |     When user right-clicks on the "statsff" area of box plot viewer
```

A misspelled menu item gets `visible menu items in context menu: General | Reset View | Markers | …`;
a phrase inside a menu that is not open gets `context menu: not open`; a misspelled property gets
the nearest captions, then all of them. A journey lists every failed scenario in that shape. There
is no matcher diff, no in-page stack and no harness line: a stack trace appears only for a
programming error in a binding, and then only its own frames. `bdd/test-results/` still has the
screenshot and the trace.

## Element phrases

A phrase resolves, in this order, at every level:

| Phrase                                     | Resolution                                                       |
|--------------------------------------------|------------------------------------------------------------------|
| `results`, `browse tab`                    | a registered element (or alias): the whole phrase wins over everything below |
| `second item …`, `last row …`, `3rd input` | an ordinal among the matches                                     |
| `save button in toolbar`                   | composition: `X in|inside|within|on|under Y` — X resolved inside Y (recursively) |
| `label of name input`, `viewers section of toolbox` | `of` names a *part* — of a registered element, or of every element of a kind (inputs: label, editor, options, error; dialogs: title, close button, footer) |
| `sequence column input`                    | a generic **kind** by its longest suffix, qualified by the rest  |
| `"Run MSA" button`, `"First name" input`   | a quoted qualifier: scope words inside it are kept, and a leading "first"/"last" is not read as an ordinal |

**Kinds** cover the whole u2 library (every `data-u2` value it stamps) plus the Dart shell's
conventions: inputs (`input`, `text input`, `text area`, `choice input`, `multi choice input`,
`number input`, `checkbox`, `date input`, `color input`, `font input`, `icon input`, `image input`,
`list input`, `map input`, `message input`, `radio input`, `slider`, `range slider`, `slider handle`,
`suggest input`, `combobox`, `multi select`, `tags input`, `file input`, `columns input`,
`column input`, `function input`, `dynamic input`, `rsa input`, `dart input`), actions (`button`,
`icon button`, `dropdown button`, `button group`, `row actions`), forms (`form`, `function form`,
`object form`, `property grid`, `property`, `category`), collections (`list`, `item`, `tree`,
`tree node`, `table`, `table row`, `virtual grid`, `functions browser`, `history browser`), display
(`icon`, `badge`, `tag`, `notification`, `progress bar`, `stat card`, `tooltip`, `tour`, `async view`,
`tray`, `heading`, `text`, `link`), navigation (`toolbar`, `menu`, `menu bar`, `menu item`,
`breadcrumbs`, `breadcrumb`), containers (`dialog`, `tabs`, `tab`, `tab panel`, `section`,
`accordion`, `accordion header`, `card`, `wizard`, `wizard step`, `splitter`, `splitter panel`,
`sash`), entities (`chip`, `entity card`, `palette`, `designer`) and the shell (`viewer`, `view`,
`element`). Each kind knows how its qualifier narrows the candidates — `data-u2-name`, the label or
title part, the text, an aria label, a placeholder, the Dart client's `name=` conventions
(`icon-scatter-plot`, `viewer-Grid`, `div-section--Viewers`, `input-host-Caption`). A row goes by
its primary text (`Abs item` finds the functions-browser row "Abs — (x) : num", `Caffeine table row`
the row whose first cell says Caffeine). Every suffix split is kept: `sequence column input` is
tried as `column input` qualified "sequence" and as `input` qualified "sequence column". `grok-bdd
lint` prints how each phrase resolved.

**Look before you phrase.** A phrase is only as good as the markup it was written against: open the
page in devtools and read the u2 contract — `data-u2` (the kind), `data-u2-name` (a deliberate
name), `data-u2-part` (label/editor/options/error), the ARIA roles and states — and write the
phrase that matches what is there (`Caffeine table row`, `Advanced pane`, `"Search users…" suggest
input`). Two traps: a qualifier that starts with an ordinal word needs quotes (`"First name"
input`), and a label equal to a kind name (`Columns input`) resolves to the labelled input first,
then to the kind. `grok-bdd lint` shows every reading the compiler kept.

**Overrides.** Register the whole phrase and composition is bypassed:
`workbench.element('save button inside toolbar', {selector: '[data-u2-name="toolbarSave"]'})`.

**Reserved names and contexts.** Element names are global and unique, and the platform base owns
the shell's: `toolbox`, `toolbox tab`, `browse tab`, `context panel`, `console`, `status bar`,
`open tableview`, `grid`. Registering one of them again is an error, even inside an app: "toolbox"
always means the Datagrok toolbox, and "grid" the grid viewer (a u2 VirtualGrid is a `virtual
grid`). A package's own names live on a *context*, a named region:

```ts
// bdd/bindings/elements.ts
import {context} from '@datagrok-libraries/bdd';
export const workbench = context('MSA workbench', {selector: '[data-u2-name="msaWorkbench"]'});
workbench.element('results', {selector: '[data-u2-name="results"]'});
workbench.element('log', {selector: '[data-u2-name="log"]', in: 'results'});

// bdd/bindings/steps.ts
import {Given} from '@datagrok-libraries/bdd';
export const openWorkbench = Given('user opens the MSA workbench', async (page) => { … },
  {enters: 'MSA workbench'});
```

After a step declared with `enters`, the context's names apply, and both they and the generic
kinds are looked up inside the context root first, then on the whole page (dialogs and
notifications are portaled out of it). Platform names keep their platform meaning; a context is an
element itself (`MSA workbench should be visible`); `leaves: true` returns to the platform
vocabulary. Choose the root with the lookup in mind: U2Demo's `U2 Demo` context is the sub-demo's
content pane, so `Tables tree node` is the page's tree, not the navigation's — which the context
names `demo navigation` for the times a feature wants it.

**Package kinds.** A package can register a generic kind of its own with `kind()` — the same
mechanism the library's vocabulary uses — for a repeated structure that carries no `data-u2`.
U2Demo's pages print their signals as "name = value" lines:

```ts
kind('readout', {
  selector: '.u2demo-status',
  match: ['label'],
  labelSelector: 'span:first-child',
  parts: {label: 'span:first-child', value: 'span:last-child'},
});
// Then value of "dose * replicates" readout should have text "750"
```

Selectors inside `labelSelector`, `parts` and `has`-style filters are evaluated *from the element*
(Playwright scopes them), so `.u2demo-status > span` would not match its own first span there;
`span:first-child` does.

## Steps

`grok-bdd list-steps` prints them all. The base vocabulary (`bindings/common/steps.ts`):

```
When  user clicks (on ){element}          user double-clicks (on ){element}     user right-clicks (on ){element}
      user hovers (over ){element}        user focuses (on ){element}           user types {string} in(to) {element}
      user enters {string} in(to) {element}   user clears {element}             user presses {key}
      user presses {key} in {element}     user selects {string} in {element}    user checks/unchecks/toggles {element}
      user opens {element}                user closes {element}                 user expands/collapses {element}
      user drags {element} to {element}   user scrolls to {element}             user navigates to {string}
      user reloads the page               user fills in:  | element | value |
Then  {element} should be/become {state}  {element} should not be/become {state}
      {element} should contain (the )text {string}      {element} should have (the )text {string}
      {element} should have (the )value {string}        {element} should have {int} item(s)/row(s)/tab(s)
      the following elements should be {state}:  | element |
Given user is logged in                   user opens {dataset} dataset          user waits for {element}
      user waits for {int} millisecond(s)   (a probe while writing a feature; a committed one needing it is missing a signal)
When  user switches to (the ){string} table view   user closes all views
      user saves the current view as project {string}   (deleted again when the feature ends)
      user opens the {string} project
```

The current table's rows and columns (`bindings/platform/data.ts`), through the JS API the way a
viewer sees them; a step that changes rows or colors snapshots every open viewer first, so the
viewer checks that follow compare with the state before it:

```
When  user clears the row selection       user deletes the selected rows
      user filters rows where {string} is between {float} and {float}    user filters out rows where {string} is {string}
      user resets the filter              user adds a calculated column {string} with formula {string}
      user removes {string} column        user removes the coloring of {string} column
      user colors {string} column linearly from {string} to {string} (over {float} to {float})
      user colors {string} column conditionally:  | range | color |     user colors {string} column categorically:  | category | color |
Then  no/some rows should be selected     all/only/some/no rows where {string} is {string} should be selected
      the table should have a current row   the table should have {int} row(s)   the table should have no rows where {string} is {string}
      {int} row(s) should pass the filter   fewer than {int} rows should pass the filter   all rows should pass the filter
      the filter should pass exactly the rows where {string} is between {float} and {float}
      table {string} should be open       table {string} should have columns {string}    table {string} should have {int} row(s)
      table {string} should have no missing values in {string} column
```

States: visible, hidden, present, absent, enabled, disabled, checked, unchecked, selected, empty,
expanded, collapsed, focused. `selected` reads whatever the element uses to say so —
`aria-selected` on options and tabs, `aria-pressed` on toggles and selectable cards,
`aria-checked` on radio-like buttons, `aria-current` on wizard steps; `expanded`/`collapsed` read
the element or the header/trigger inside it (a section, an accordion pane, a property category); a
phrase matching several elements (stacked balloons) is `visible` when any is and `hidden` when none
is. Outcomes are Playwright's retrying `expect`, so a feature never needs a wait step. When two
definitions match a step, the one that spells out more of the text wins, then the one with fewer
parameters; an exact tie is a compile error, never a silent pick. A gesture on a phrase that
matches several elements acts on the visible ones (a Dart menu keeps hidden mirrors of its items);
several visible matches are reported, never picked.

`user selects {string} in {element}` takes a native `<select>` as it is; anything else gets its
editor clicked (a combobox or typeahead also gets ArrowDown, since those open on a keystroke) and
the option picked from the popup by its whole text, its primary-text part, its title or aria label
(icon cells), then by substring. `user expands/collapses` clicks the element's own `aria-expanded`
control, or a tree row's twistie. `user types`/`enters` and the value checks address the editor:
the input inside an input root, the contenteditable of a message input, the trigger of a picker,
the element itself when nothing is inside (a slider handle, a list).

A step definition is an exported `const` (the generated spec imports it by name):

```ts
export const openDataset = Given('user opens {dataset} dataset', async (page: Page, dataset: DatasetEntry) => {
  await openTableFromFile(page, dataset.path);      // from @datagrok-libraries/test
}, {tier: 'api', description: 'OpenFile through the JS API — provenance as in the UI'});
```

(The real one also waits for the platform's semantic-type-detected event for that table: a step is
over when the platform is done with it, so background work never lands on the next step.)

Parameter types: `{element}` (any phrase), `{widget}` (a phrase naming a viewer or a widget — the
`viewers` tier takes nothing else where it reads or changes one), `{dataset}` (a registered alias or a
platform path), `{viewer}` (a viewer type), `{key}` (`Enter`, `Control+A`, `ArrowDown`), `{state}`, plus Cucumber's `{string}`,
`{int}`, `{float}`. Runtime helpers for your own steps come from `@datagrok-libraries/bdd/runtime`:
`locate(page, el('…'))`, `gestures.*`, `expectState`, `expectText`, `expectValue`, `expectCount`.

## Tiers

`bindings/common` and `bindings/platform` load for every project. Vocabulary that only some
packages need is a tier under `bindings/tiers/<name>/`; a package opts in with
`{"tiers": ["viewers"]}` in `bdd/bdd.config.json`. A new tier is a directory of binding modules
like any other; nothing to register.

### The `viewers` tier

Viewers on the current table view, written the way the platform sees them — properties by their
caption, context menus by their path, canvas regions by the names the viewer reports, repaints by
the viewer's own render event. The first features written with it are the six box plot journeys
under `packages/UsageAnalysis/bdd/features/viewers/` (property surface, group comparison,
selection, filter, statistics and coloring, settings ladder): 50 scenarios, the whole of six
hand-written Playwright specs and their helpers, in 63 s.

```gherkin
Background:
  Given user is logged in
  And user opens demog-1000 dataset
  And user adds a box plot viewer with:
    | Value      | AGE |
    | Category 1 | SEX |

Scenario: Context menus as property paths
  When user picks "Misc > Show Inside Values" from the context menu of box plot viewer
  Then "Show Inside Values" property of box plot viewer should be "false"
  And box plot viewer should have less ink than before
  When user sets "Marker Size Column" property of box plot viewer to "WEIGHT"
  And user opens the context menu of box plot viewer
  And user hovers over Markers menu item in context menu
  Then "Markers > Size" menu item in context menu should be disabled
```

```
Given user adds (a ){viewer} viewer                 user adds (a ){viewer} viewer with:  | caption | value |
      user listens for {string} event on {widget}  user switches to (the ){string} table view
When  user sets {string} property of {widget} to {string}       user sets properties of {widget}:  | caption | value |
      user picks {string} from the context menu of {element}
      user picks {string} from the context menu of the {string} area of {widget}
      user opens the context menu of {element}      user right-clicks on the {string} area of {widget}
      user closes the context menu                  user clicks / double-clicks / hovers over the {string} area of {widget}
      user clicks on the {string} area of {widget} holding {key}    user drags a selection box over the {string} area of {widget}
      user moves the pointer away from {element}    user resizes {widget} to {int} by {int}    user resizes {widget} to {int} wide
      user restores the size of {widget}            user takes a snapshot of {widget}          user remembers the value range of {widget}
      user saves the layout of the current table view   user loads the saved layout
Then  {string} property of {widget} should be {string}          {string} property of {widget} should not be {string}
      properties of {widget} should be:  | caption | value |
      {widget} should have repainted                {widget} should not have repainted           {widget} should be painted
      {widget} should have less/more ink than before    {widget} should be painted in at least {int} colors
      the {string} area of {widget} should be painted   {widget} should (not )have a(n) {string} area
      {widget} should show more/less selection highlight than before    {widget} should show a/no selection highlight
      {widget} should show a narrower/wider value range than before     {widget} should show the same value range as before
      the value range of {widget} should lie within {string} column     {widget} should show the remembered value range
      {string} event should have fired on {widget}  {string} event should not have fired on {widget}
      no errors should have been logged             no error or warning balloon should have been shown
      the tooltip should show columns {string}      the tooltip should not show columns {string}
```

A property is named by its caption as the property panel shows it (`Value`, `Category 1`,
`Show Markers`, `Marker Size Column`) or by its name (`showInsideValues`); values are `true`/`false`,
numbers, `#rrggbb` for colors (read back the same way), `""` for none, `\n` for a line break. A
menu path is `"Group > Item"`. A hit area is a name the viewer reports (`grok-bdd lint` cannot
list them yet; the box plot has `view`, `x axis`, `y axis`, `stats`, `p value`, `group comparison`,
`color scale`, `marker`, `category <label>` and `<label> values` per category, `p value of <group>`
and `<effect> effect` under group comparison; the bar chart `view`, `x axis`, `y axis`, `bar
<category>`). Every property set, menu pick, area click, hover and resize snapshots the canvas,
the selection-colored pixels and the value range first, so `should have repainted`, `less/more
ink`, `more/less selection highlight` and `a narrower/wider/the same value range than before`
compare with the state before the last change; the data steps snapshot every viewer. `no errors
should have been logged` is the page's console errors and uncaught exceptions since the previous
check (or the login step); a resource the stand does not serve is not an error. `no error or
warning balloon should have been shown` reads the platform's `d4-balloon-shown` events the same way.

Viewers on a bdd page render immediately — `viewer.immediateRendering` is set on every viewer the
page holds or adds — so nothing in the tier sleeps: a change is followed by the viewer's render
event, a context menu by `onContextMenuShown`. **When a step would need a wait, the platform is
missing a signal; it goes into the core, not into the step** (the library's `CLAUDE.md` keeps the
list of what was added that way). `user listens for {string} event on {widget}` subscribes to the
viewer's event once; `should have fired` reads the count and ends the subscription, and a viewer
that closes drops its subscriptions itself.

## Generated specs

One `test()` per scenario and per outline row (`Every method reports itself [method=muscle]`), one
`test.step` per Gherkin step, tags as Playwright tags, `@realizes:<feature>` tags collected into a
`sub_features_covered` header. A `@journey` feature is one `test()` with the Background inline and
every scenario a `run.scenario(...)` soft step, closed by `run.finish()`. Library bindings are imported by package subpath, the package's own
by relative path, and phrases are emitted as names (`el('…')`, `ds('…')`), never selectors — so a
selector fix never regenerates anything, and `grok-bdd compile --check` fails only when a feature
changed without its spec.

## VS Code

The official Cucumber extension (`cucumberopen.cucumber-official`) needs, in the package's
`.vscode/settings.json`: `cucumber.features` → `bdd/features/**/*.feature`, `cucumber.glue` →
`bdd/bindings/**/*.ts` and `node_modules/@datagrok-libraries/bdd/bindings/**/*.ts`, and one entry in
`cucumber.parameterTypes` for `state` (the extension reads the other custom types — element, dataset,
viewer, key — from the `defineParameterType` calls in the glue itself; listing them twice only
produces "already a parameter type" errors). See `packages/U2Demo/.vscode/settings.json`. A step
shown as undefined while `grok-bdd lint` resolves it means the extension's view differs from the
compiler's: check the settings file is valid JSON (a `\\w` that became `\w` silently turns
`[\w -]+?` into `[w -]+?`) and that the glue globs reach the tier directories.

## Developing the library

`npm run build` compiles `src/`, `bindings/` and the Playwright config to `dist/` (what the
`exports` map and the command use); `npm run test:unit` runs the engine tests; the library is a
project itself (`features/platform`, tier `viewers`) — `npm test` builds and runs it.

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
  first pointer move, and the semantic-type detection changes; the package features (Bio) on
  2026-09-08's `Func.topMenu`, `aria-disabled` on dialog buttons and the grid's
  `getWidgetStatus`. On an older platform the `viewers` features fail on those steps; the u2
  features do not care.
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
grok-bdd run [--headed] [-g "name"] [--reporter=list]   # compile --check, then Playwright (headed runs the same: the browser launches with its accelerated 2D canvas off, so a repaint check reads the same pixels as headless)
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

One browser page per worker: the first scenario the worker runs opens it and boots the shell
(about 4 s), every scenario ends with the shell reset (dialogs and popups closed,
`grok.shell.closeAll()`, the Home view current), and every later feature starts on that reset
shell — `user is logged in` only navigates when the page is not in the shell yet, and a
package's init step is free once the package is up. What a feature leaves on the server it puts
back itself (`atFeatureEnd`). Folders group the features by subject (`features/viewers/box-plot/`)
and are the unit of `grok-bdd run generated/<folder>`. Playwright still runs and reports one test
per scenario (and per outline row), so `-g`, tags, retries and traces work as usual, each test
with its own trace; a trace keeps the actions, console and network, and a failed test its
screenshot — `grok-bdd run --trace on` records DOM snapshots and a screenshot per action too,
`--video on` a video. A `Background` runs before every scenario, as Gherkin says.

**`@journey`** on the feature changes that: the feature is one test, the Background runs once, and
the scenarios run in order on the same shell state, each a soft step — a failing scenario is
recorded and the next one still runs, and the test fails at the end listing them. Use it for a
property surface walked section by section, where re-opening the data and the viewer thirteen
times would cost more than the checks (the box plot: 13 scenarios in 14 s, 6 of them the login and
the open). Each scenario then puts back what it changed, so the next starts where the Background
left off, and owns its error and balloon floors — what an earlier scenario logged is not charged to
it. `-g` selects the whole journey; the report shows every scenario and step under it.

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
      user clicks on {element} holding {key}   (Control adds to a selection, Shift extends it, Control+Shift removes)
      user pastes {string} in(to) {element}    (through the clipboard and Ctrl+V, so the platform's paste handling runs)
Then  {element} should be/become {state}  {element} should not be/become {state}
      {element} should contain (the )text {string}      {element} should have (the )text {string}
      {element} should have (the )value {string}        {element} should have {int} item(s)/row(s)/tab(s)
      the following elements should be {state}:  | element |
Given user is logged in                   user opens {dataset} dataset          user waits for {element}
      user waits for {int} millisecond(s)   (a probe while writing a feature; a committed one needing it is missing a signal)
When  user switches to (the ){string} table view   user switches to (the ){string} view   user closes all views
      user saves the current view as project {string}   (deleted again when the feature ends)
      user opens the {string} project
```

The top menu and what its commands do (`bindings/platform/commands.ts`): a path is picked by
real pointer moves over the Dart menu bar (`"Bio > Analyze > MSA..."`, the labels as shown), the
function call the command starts is watched — the platform announces every call, and the one
registered under the picked path (`Func.topMenu`) is the command — and the columns it adds are read
against the columns the table had when it was picked:

```
When  user picks {string} from the top menu     user opens {string} in the top menu    user closes the top menu
Then  the top menu should list:                 (a table of paths; each group opens once for all its leaves)
      the top menu command should have completed          (its function call has ended, dialog and all; up to two minutes)
      {int} new column(s) should have been added          a new column {string} should have been added
      a new column matching {string} should have been added   no new column should have been added
```

The columns of the current table as facts (`bindings/platform/columns.ts`), a package function
and its result (`bindings/platform/functions.ts`), and a dataset opened by its first rows:

```
Given user opens {dataset} dataset keeping the first {int} rows      … keeping the first {int} rows as {string}
Then  {string} column should have semantic type {string}    {string} column should have units {string}
      {string} column should have tag {string} equal to {string}    {string} column should have type {string}
      {string} column should have no/missing values     every value of {string} column should match/contain {string}
      every value of {string} column should lie between {float} and {float}
      every value of {string} column should have the same length (within each {string} value)
      every value of {string} column should be {string} and {string} of the same row joined by {string}
      the value of {string} column in row {int} should be {string}    {string} column should have its maximum in row {int}
      {string} column should have at least {int} distinct values
      the table should have a column {string}   the table should not have a column {string}   the table should have {int} column(s)
When  user makes row {int} current     user makes the last row current     Then  row {int} should be current
When  user calls {string} function     user calls {string} function with:  | name | value |   (column:X, table, numbers, true/false)
Then  the result should be empty / contain text {string} / match {string} / have methods {string} / have a {string} of {string}
      the result should be a list of {int} or more items    the result should be a table with columns {string}
      every column of the result table should be filled in row {int}
Then  the filter panel should have {int} filter(s)    the filter panel should have a filter on {string} column
      the filter should pass exactly the rows where {string} contains {string}
      only rows where {string} starts with {string} should be selected
      the {string} view should be current     When user closes the current view
```

A package's word about work finished off-screen, a file chooser, the clipboard, an app, and the
grid (`bindings/platform/events.ts`, `common/steps.ts`, `platform/steps.ts`, the viewers tier):

```
Given user listens for {string} custom event        (grok.events.onCustomEvent — Bio's bio-monomer-lib-loaded)
Then  the {string} custom event should have fired   the {string} custom event should not have fired
When  user uploads {string} through {element}       (a file of the bdd project, "fixtures/lib.json", into the chooser the element opens)
Then  the clipboard should contain/have (the )text {string}
Given user opens the {string} app                   (the app function by its name; done when its view is current)
Then  the {string} area of {widget} should be painted in at least {int} colors   (hues, greys and white aside)
```

`grid` is a widget: the platform's grid reports every visible cell as `cell <row> of <column>`
(rows as the table counts them, from 1), `row header <row>` and `header <column>` hit areas, and
the readings `rows shown`, `rows`, `columns shown`, `cell type of <column>` for each visible
column, `current row` and `current column` — so `the "cell type of fasta" reading of grid should
be "sequence"`, `user picks "Copy > helm" from the context menu of the "cell 1 of fasta" area of
grid`, `user clicks on the "cell 2 of fasta" area of grid`, and
`the "cell 1 of HELM string" area of grid should be painted in at least 3 colors` for a renderer
that colors its monomers.

The current table's rows and columns (`bindings/platform/data.ts`), through the JS API the way a
viewer sees them; a step that changes rows or colors snapshots every open viewer first, so the
viewer checks that follow compare with the state before it:

```
When  user clears the row selection       user deletes the selected rows
      user selects rows where {string} is {string}     user selects rows where {string} is one of {string}
      user filters rows where {string} is between {float} and {float}    user filters out rows where {string} is {string}
      user filters rows where {string} is {string}     user resets the filter
      user adds a calculated column {string} with formula {string}
      user removes {string} column        user removes the coloring of {string} column
      user colors {string} column linearly from {string} to {string} (over {float} to {float})
      user colors {string} column conditionally:  | range | color |     user colors {string} column categorically:  | category | color |
      user colors {string} column {word} again   (switches the type back on, keeping the colors it stored)
      user colors {string} column linked to {string} column    user colors the text of {string} column linked to {string} column
      user applies the coloring of {string} column to {string} column    user inverts the color scheme of {string} column
      user sets {string} column in row {int} to {string}   ("NaN", "Infinity", "-Infinity" and "" on a numeric column)
      user selects rows where {string} is between {float} and {float}
Then  no/some rows should be selected     all/only/some/no rows where {string} is {string} should be selected
      all/no rows where {string} is {string} should pass the filter
      the table should have a current row   the table should have {int} row(s)   the table should have no rows where {string} is {string}
      {int} row(s) should pass the filter   fewer than {int} rows should pass the filter   all rows should pass the filter
      the filter should pass exactly the rows where {string} is between {float} and {float}
      the filter should pass exactly the rows where {string} is {string}
      table {string} should be open       table {string} should have columns {string}    table {string} should have {int} row(s)
      table {string} should have no missing values in {string} column
      {string} column should have no color coding    {string} column should be color-coded categorically
      {string} column should be color-coded linearly/conditionally/categorically/linked
      the coloring of {string} column should be linked to {string} column    the text of {string} column should be color-coded
      the color scheme of {string} column should be {string}   (the stops of a linear scheme in order)
      {int} row(s) should be selected     rows {int} to {int} should be selected     all rows should be selected
      columns {string} should be selected    no columns should be selected    every selected row should pass the filter
      {int} rows of table {string} should pass the filter    {int} rows of table {string} should be selected
      {string} column should be true exactly where the filter passes
      the categorical color of {string} in {string} column should be {string}
When  user selects the first {int} rows
      user links table {string} to table {string} as {string}:  | key in the first | key in the second |   ("filter to filter", "selection to filter", …)
```

The filter panel through its own API (`bindings/platform/data.ts` too): the panel is opened, a
card's criterion set as a category click or a range drag would leave it, and its state read back —
the gestures on the cards are the filter panel features' own subject, these are their setup:

```
When  user opens the filter panel      user opens an empty filter panel
      user adds a categorical filter on {string} keeping {string}     user adds a range filter on {string} from {float} to {float}
      user configures the hierarchical filter with columns {string}
Then  the filter on {string} column should keep only {string}
      the filter on {string} column should be filtering     the filter on {string} column should not be filtering
```

`filter panel` names the Filters viewer (`[name="viewer-Filters"]`) with the parts `counter`,
`master`, `search`, `add filter selector`, `reset icon`, `search icon` and `expand icon`
(`counter of filter panel should be hidden`), and `{widget}` accepts it, so its readings and hit
areas read as a viewer's; a card is a `filter card` by its caption (`"RACE" filter card`, parts
`caption`, `checkbox`, `indicator`, `mode`, `summary`, `body`, `close`, `search icon`, `search`),
`disabled` while suspended. A hierarchical node is `partially checked` when its children
disagree; an input the platform refuses is `invalid`.

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

(The real one also waits for the platform's semantic-type-detected event for that table, and makes
row 0 current itself — the view would do that a second after the grid appears, repainting every
viewer on whatever step is running then: a step is over when the platform is done with it, so
background work never lands on the next step.)

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
under `packages/UsageAnalysis/bdd/features/viewers/box-plot/` (property surface, group comparison,
selection, filter, statistics and coloring, settings ladder): 50 scenarios, the whole of six
hand-written Playwright specs and their helpers, in 36 s on one page, backward-matched against
the specs they replaced (see "Translating an existing spec" below).

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
Then  the open menu should (not )list {string}      ("Group > Item" — the groups are opened to reach the item)
      user clicks on the {string} area of {widget} holding {key}    user drags a selection box over the {string} area of {widget}
      user drags a selection box from the {string} area to the {string} area of {widget}
      user drags across the {string} area of {widget}    user scrolls the mouse wheel up/down over the {string} area of {widget}
      user drags the {string} area of {widget} to the {string} area    user drags the {string} area of {widget} by {int} pixels to the left/right/up/down
      user drags a deselection box over the {string} area of {widget} (Control+Shift)    user drags a zoom box over the {string} area of {widget} (Alt)
      user enters {string} into the {string} area of {widget}    (an area that holds an editor: click, select all, type, Enter)
      user remembers the {string} reading of {widget}    user clicks on {string} item in the legend of {widget} (holding {key})
      user clicks on the cross of {string} item in the legend of {widget}    user drags the legend splitter of {widget} by {int} pixels
      user moves the pointer away from {element}    user resizes {widget} to {int} by {int}    user resizes {widget} to {int} wide
      user restores the size of {widget}            user takes a snapshot of {widget}          user remembers the value range of {widget}
      user saves the layout of the current table view   user saves the layout of the current table view to the server   user loads the saved layout
Then  {string} property of {widget} should be {string}          {string} property of {widget} should not be {string}
      properties of {widget} should be:  | caption | value |    {widget} should be bound to table {string}
      {widget} should have repainted                {widget} should have repainted by at least {int} pixels
      {widget} should not have repainted            {widget} should be painted
      {widget} should have less/more ink than before    {widget} should be painted in at least {int} colors
      the {string} area of {widget} should be painted   the {string} area of {widget} should have less/more ink than before
      the {string} area of {widget} should have repainted
      the open tableview should have {int} {viewer} viewer(s)
      the {string} area of {widget} should contain the color {string}
      the {string} and {string} areas of {widget} should be painted in different colors
      the {string} and {string} areas of {widget} should be painted in the same colors
      the {string} area of {widget} should not contain the color {string}
      the {string} area of {widget} should be at least {int} pixels tall/wide
      the {string} area of {widget} should be taller/wider/shorter/narrower than before
      {widget} should (not )have a(n) {string} area
      the {string} reading of {widget} should (not )be as remembered
      the {string} reading of {widget} should not be {string}    the {string} reading of {widget} should be a finite number
      the {string} and {string} readings of {widget} should be the same / should differ
      the legend of {widget} should list {int} item(s) / fewer items than before / the same items as before
      the legend of {widget} should be docked / in a corner / collapsed to the mini icon / shown in the tooltip
      the legend of {widget} should (not )be in the {string} slot    the legend of {widget} should be placed as before
      the {string} item in the legend of {widget} should be colored {string}
      the {string} and {string} items in the legend of {widget} should be colored differently
      the tooltip should show {string} as {string}
      {widget} should show more/less selection highlight than before    {widget} should show a/no selection highlight
      {widget} should show a narrower/wider value range than before     {widget} should show the same value range as before
      the value range of {widget} should lie within {string} column     {widget} should show the remembered value range
      the color scale of {widget} should cover a narrower/wider range than before
      the color scale of {widget} should cover the same range as before
      the {string} reading of {widget} should be {float}    the {string} reading of {widget} should be lower/higher than before
      the {string} reading of {widget} should differ from before    the {string} reading of {widget} should be the same as before
      legend of {widget} should be visible/hidden / have {int} items / contain text {string}    the legend of {widget} should be on the left/right/top/bottom
      {string} event should have fired on {widget}  {string} event should not have fired on {widget}
      no errors should have been logged             no error or warning balloon should have been shown
      an error/a warning balloon should have been shown     an error/a warning balloon containing {string} should have been shown
      {string} property of {widget} should contain {string}
      the {string} reading of {widget} should be {string}   the {string} reading of {widget} should be at least {float}
      the tooltip should show columns {string}      the tooltip should not show columns {string}
      the tooltip should show some columns
```

A property is named by its caption as the property panel shows it (`Value`, `Category 1`,
`Show Markers`, `Marker Size Column`) or by its name (`showInsideValues`); values are `true`/`false`,
numbers, `#rrggbb` for colors (read back the same way), `""` for none, `\n` for a line break. A
menu path is `"Group > Item"`. A hit area is a name the viewer reports (`grok-bdd lint` cannot
list them yet; the box plot has `view`, `x axis`, `y axis`, `stats`, `p value`, `group comparison`,
`color scale`, `marker`, `category <label>` and `<label> values` per category, `p value of <group>`
and `<effect> effect` under group comparison; the bar chart `view`, `x axis`, `y axis`, `bar
<category>` and `bar <category> | <stack>`; the 3D scatter plot `view` and `point`; under group
comparison with a control, `control band` or `control band <stratum>`). A reading is a name the
viewer reports in `getWidgetStatus().values`: every viewer with a row filter `rows shown`, the box
plot `color scale min` / `max`, the bar chart `bars`, `stack segments`, `clipped bars`, the 3D
scatter plot `camera x` / `y` / `z` / `camera distance` and `scene signature` — its WebGL canvas
has no pixels to read, so `the "scene signature" reading … should differ from before` is its
`repainted`. `legend of {widget}` is the viewer's legend element (`[name="legend"]`, its rows the
legend items); `should have {int} items` counts the rendered rows (the legend virtualises long
lists), `should list {int} items` the total the legend publishes, and its mode, slot and keys come
from the same `data-legend-*` attributes the legend writes on every commit — the library adds
nothing to the legend's render path. A `legend item` is a kind (`"R_ONE" legend item in legend of
scatter plot viewer`, `selected` while its category filters the viewer, parts `label`, `cross`,
`thumbnail`, `marker`); the legend's chrome are viewer parts (`mini legend icon`, `legend
splitter`, `legend inner splitter`, `legend close chevron`, `legend markers selector`), the hover
pickers page elements (`color picker icon`, `marker picker icon`). A viewer with no canvas (a
form, the filter panel) reports its hit areas relative to its root and has no pixel steps; a
viewer with an `overlay` canvas of the same size as its canvas (the scatter plot's regression lines
and labels, the grid's selection) has both composited into every pixel reading.
Every property set, menu pick, area click, hover and resize snapshots the canvas, the ink of every
hit area, the selection-colored pixels, the value range and the color scale's range first, so
`should have repainted`, `less/more ink` (of the canvas or of one area), `more/less selection
highlight`, `a narrower/wider/the same value range than before` and `the color scale … than before`
compare with the state before the last change; the data steps snapshot every viewer, and return
once every viewer has drawn the change.

Say what the claim is. `should have repainted` is a change detector — any pixel — and proves a
toggle reached the canvas, nothing about what was drawn; a chrome toggle with no shape of its own
(an axis, a selector) takes `by at least N pixels`. A shape gets its own evidence: `the "M values"
area … more ink than before` for a violin, `should have a "stats" area` for the strip, `should
contain the color` for a coloring, `the "M values" and "F values" areas … different colors` for
per-category hues, `should be bound to table` for a rebind (the Table property is what was asked,
`viewer.dataFrame` what happened). `more/less selection highlight` needs a margin the selection
warrants (a floor of 200 device pixels, two per selected row, capped at a quarter of the view) —
one more orange pixel is the hover halo, not a selection. The negative checks — `should not have
repainted`, `the same value range as before` — read once the viewer is quiet, so a reset that
lands a tick after the change fails them instead of slipping past. `no errors should have been
logged` is the page's console errors and uncaught exceptions since the previous check, the
scenario's start (a journey scenario owns its floor) or the login; a resource the stand does not
serve is not an error. `no error or warning balloon should have been shown` reads the platform's
`d4-balloon-shown` events the same way.

A JS viewer takes part on the same terms by giving the runtime what a Dart viewer gives it:
`getWidgetStatus()` with its canvas under `parts`, `hitAreas` in CSS px of that canvas and named
`values`; a `get isRenderPending()` that is true from a render request to the paint; and an
`onRendered` observable that fires after every render pass (the host's `onViewerRendered` never
fires for a JS viewer). Bio's WebLogo (`position <label>`, `monomer <M> at position <label>`,
`positions shown`, `rows shown`, `rows selected`) and its similarity and diversity search viewers
(`target row`, `neighbours`, `neighbour set`; `subset size`, `subset`, `distinct sequences`;
both `source column` and `limit`) are the first (`packages/Bio/src/viewers/web-logo-viewer.ts`,
`src/analysis/sequence-search-base-viewer.ts`).

Viewers on a bdd page render immediately — `viewer.immediateRendering` is set on every viewer the
page holds or adds — so nothing in the tier sleeps: a change is followed by the viewer's own word
that nothing is pending any more (`viewer.isRenderPending`: a debounced refresh armed or a repaint
requested), a context menu by `onContextMenuShown`. A property that paints nothing costs no wait,
a repaint is waited for as long as it takes, and a repaint that never lands is reported as the
platform failure it is. **When a step would need a wait, the platform is missing a signal; it goes
into the core, not into the step** (the library's `CLAUDE.md` keeps the list of what was added that
way). `user listens for {string} event on {widget}` subscribes to the viewer's event once; `should
have fired` reads the count and ends the subscription, and a viewer that closes drops its
subscriptions itself.

**Translating an existing spec** into a feature, and proving the feature tests what it claims, is
the `/bdd-translate` skill (`public/.claude/skills/bdd-translate/SKILL.md`): translate, then have
one independent reviewer per old-spec/feature pair backward-match every old assertion and try to
break every new one, then fix in the core, the library and the feature, in that order. The first
round on the box plot found the checks that were green for the wrong reason (a one-pixel
"repainted", a tooltip that kept its last text, a "same range" read before the reset landed, a
settle cap that hid a late repaint) and fixed each where it belonged.

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

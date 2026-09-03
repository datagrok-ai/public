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
npm i -g @datagrok-libraries/bdd   # the `grok-bdd` command (until it is on npm: `npm link` from libraries/bdd)
cd <package>
grok-bdd init                      # bootstraps bdd/ — see below — and the manifest entries
npm i                              # installs @datagrok-libraries/bdd and @playwright/test as dev dependencies
grok-bdd run                       # the smoke feature: logged in, Browse visible
```

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
grok-bdd compile            # features/** → generated/**  (+ diagnostics)
grok-bdd compile --check    # fails when a committed spec is stale — the CI gate
grok-bdd lint               # diagnostics only
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
```

States: visible, hidden, present, absent, enabled, disabled, checked, unchecked, selected, empty,
expanded, collapsed, focused. `selected` reads whatever the element uses to say so —
`aria-selected` on options and tabs, `aria-pressed` on toggles and selectable cards,
`aria-checked` on radio-like buttons, `aria-current` on wizard steps; `expanded`/`collapsed` read
the element or the header/trigger inside it (a section, an accordion pane, a property category); a
phrase matching several elements (stacked balloons) is `visible` when any is and `hidden` when none
is. Outcomes are Playwright's retrying `expect`, so a feature never needs a wait step. When two
definitions match a step, the one with fewer parameters wins, then the one that spells out more of
the text; an exact tie is a compile error, never a silent pick.

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

Parameter types: `{element}` (any phrase), `{dataset}` (a registered alias or a platform path),
`{viewer}`, `{key}` (`Enter`, `Control+A`, `ArrowDown`), `{state}`, plus Cucumber's `{string}`,
`{int}`, `{float}`. Runtime helpers for your own steps come from `@datagrok-libraries/bdd/runtime`:
`locate(page, el('…'))`, `gestures.*`, `expectState`, `expectText`, `expectValue`, `expectCount`.

## Tiers

`bindings/common` and `bindings/platform` load for every project. Vocabulary that only some
packages need is a tier under `bindings/tiers/<name>/` (today: `viewers` — toolbox and viewer
steps); a package opts in with `{"tiers": ["viewers"]}` in `bdd/bdd.config.json`. A new tier is a
directory of binding modules like any other; nothing to register.

## Generated specs

One `test()` per scenario and per outline row (`Every method reports itself [method=muscle]`), one
`test.step` per Gherkin step, tags as Playwright tags, `@realizes:<feature>` tags collected into a
`sub_features_covered` header. Library bindings are imported by package subpath, the package's own
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

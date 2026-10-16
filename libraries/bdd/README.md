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

**What is not a feature.** A scenario that runs a server-side Python/R/Julia script, a Jupyter kernel
or a Docker container, or reaches an outside web service, or that has nothing UI-specific (a
function called and its result checked),
is a package test or an `ApiTests` test, never a feature. The same goes for a TestTrack case marked
`manual-only` or `apitest`. The rule and its reasons are in `CLAUDE.md`, "What never becomes a feature".

**A stand that lacks a capability skips, it does not fail.** The suites run on CI, dev, local stands
and public, which differ: `Given the stand runs the "Jupyter" service` and `Given the stand can reach
the database of the "<name>" connection` are gates placed right before the first step that needs the
capability — in the Background only when every scenario needs it from its first step — so the steps before
them run everywhere and fail as usual,
and the rest of the test is reported skipped with the reason when the stand has not got it. The service
gate goes by the health the stand reports; a stand that reports none (a dev stack) lets the test go on.
`Given the "Chem" package is installed` gates on a package the feature needs but does not test; it
may stand in the Background when every scenario needs it, and never names the package under test.
`Given the stand has a reachable "<name>" connection` is the connection gate for a connection a package
brings rather than the feature (absent skips too), and `Given the stand serves the help pages` gates a
claim on the help panel's content.

**Nothing stays on the server.** Whatever a feature adds or changes on the server — entities, files,
database rows, the layout or chat the UI makes on the side, a setting or configuration of something
it does not own — is removed or restored at feature end and swept again at its start, and the
cleanup proves it is gone. See `CLAUDE.md`, "Everything a feature puts on the server goes".

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

A slow stand may need longer budgets: `BDD_EXPECT_TIMEOUT` raises what every check waits (15 s)
and `BDD_COMMAND_TIMEOUT` how long a top-menu command has to add its columns (180 s).
`BDD_FRESH_PAGE=1` reloads the shell before every feature instead of resetting it, which tells a
feature that fails on what an earlier one left behind. `BDD_MOLECULE_SKETCHER=<name>` runs the
suite with one molecule sketcher: `Given the molecule sketcher is "…"` pins that one instead of the
feature's own, and a feature about one sketcher's own controls, tagged `@sketcher-controls`, skips
with its reason under an override that names another.

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
features rely on signals the core gained for them), and with `TableView.getSnapshot` skipping a
hidden view (without it the Save project dialog logs "Unable to find element in cloned iframe" once
per table view not shown, and a feature that saves a project with several views fails its error
check); a login — global setup mints a token from the
dev key of the `localhost` entry in `~/.grok/config.yaml` (or `DATAGROK_SERVER=<name>`), falls back
to the login form with `DATAGROK_LOGIN`/`DATAGROK_PASSWORD`, and a CI runner passes
`DATAGROK_AUTH_TOKEN`; another stand through `DATAGROK_URL=https://…`. The datasets the features
open are registered in `bindings/platform/datasets.ts`: `demog-1000` is
`System:DemoFiles/demog-1000.csv`, uploaded once with
`grok s files put public/packages/ApiTests/files/datasets/demog-1000.csv "System:DemoFiles/demog-1000.csv" --host localhost`;
`spgi` comes with the published `Chem` package (`grok s packages install Chem`). A sharing feature
shares with the account in `DATAGROK_SHARING_LOGIN`, or, unset, with the `bddsecond` user the
setup creates on the stand (a dev key is needed for that; users cannot be deleted, so it stays). A
feature that signs in as that account or another does it on its own page, and the running account
is back when the feature ends; a feature that empties the second account's notifications runs only
against a `bdd…` fixture account.

`grok-bdd init` runs in the package directory and creates what is missing, never overwriting:

```
<package>/bdd/
  package.json        {"type": "module"}   — the specs and bindings here are ES modules
  bdd.config.json     {"tiers": ["viewers"]}   — optional: library tiers beyond the base, the audits' modes
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
shell — until the page has run `BDD_PAGE_MAX_FEATURES` features (25) or its renderer holds more than
`BDD_PAGE_MAX_MB` (3000): closed views stay in the renderer's memory, so the next feature then opens a
new page in the same context. A dataset is read from the server once per page and every feature gets a clone, with the
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
`@known-failure:<line>` pins the defect to the step at that feature line: a failure of any other step
fails the test as a plain failure would. Pin every tag — a scenario that breaks earlier for another
reason stays green under a bare one.

**`@skip`** on a scenario, or on a journey's feature, declares the test and does not run it; the run
reports it skipped. It is for a test the lead has set aside (a flaky feature of an app nobody uses),
never for a defect (`@known-failure`) or a stand without a capability (the gates). The feature's
description says, in a `TODO:` paragraph, why and what brings it back.

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

Every run also leaves Playwright's JSON report in the project's `test-results/report.json`, with
the stand and the machine it ran on: the failed step of a test is its deepest step with an `error`.
A `--reporter` given on the command line gets `json` added, to the same file; one that names `json`
itself (`--reporter=list,json`) sends the report to stdout, unless `PLAYWRIGHT_JSON_OUTPUT_FILE` names
a file.
The run history of UsageAnalysis (`packages/UsageAnalysis/bdd/history`) keeps such reports, when
asked to, as dated records with a page to follow times, failures and flakes over them.

## Audits

Every scenario that passed as expected on a page showing u2 (in the current view, or a u2 dialog,
popup or menu over it) is audited by the harness, with no step in the feature:

- **layout** — the scenario's final state, before the shell reset, through the layout audit over
  the visible u2 controls (`src/runtime/layout.ts`): a control that shows something but has no size
  or is out of reach, two controls lying on one another (compared within one layer: a dialog never
  overlaps the page under it), text cut off or ellipsized, content spilling out of its control, a
  label off its editor, form editors out of their column. No screenshot is compared; a busy control
  (`aria-busy`) is waited for up to 2 s, then left out.
- **errors** — the console errors and uncaught exceptions the scenario logged, before the reset too.
- **leaks** — after the reset has closed the view: what u2 created since the scenario started and
  still holds live — scopes, timers, global listeners, each with its creation stack (`DG.U2.debug`).
  u2 tracks them only in the pages of a project that audits leaks (`__U2_DEBUG__`, set before the
  platform loads); a scenario that starts on a page not yet in the shell has no baseline and no check.

A `@journey` audits layout and errors as each scenario ends and leaks once, at the end of the test;
a test with a `@known-failure` scenario skips its end-of-test audits, leaks included. The modes —
`fail`, `report` (printed and annotated on the test as `bdd-audit`) or `off` — and the allowlists
are the project's, in `bdd.config.json`:

```json
{"tiers": ["viewers"],
 "audits": {"layout": "fail", "errors": "fail", "leaks": "report",
   "allow": {
     "layout": [{"rule": "overlap", "selector": ".my-pkg-slider", "reason": "the thumb overlaps the track on purpose"}],
     "errors": [{"message": "^Unable to get project asset", "url": "my-pkg", "reason": "the stand has no Chem"}],
     "leaks": [{"label": "interval timer", "stack": "pollQueue", "reason": "the queue poll outlives the view by design"}]}}}
```

Without the block, layout and errors are `report` and leaks `off`; `BDD_AUDIT=<mode>` sets all three
for a run, `BDD_AUDIT_LAYOUT`, `BDD_AUDIT_ERRORS` and `BDD_AUDIT_LEAKS` one each. The suite runs in a
1920×1080 window; `BDD_VIEWPORT=<w>x<h>` (say `1280x800`) runs it at another size, so the layout audit
checks every scenario's end state on a smaller screen. Every allowlist entry
says why the finding is no defect (`reason`, required): a layout entry drops the findings of its rule
whose element is inside `selector`; `message`, `url`, `label` and `stack` are regular expressions, a
leak's label being the scope's label, `<kind> timer` or `<target> <type> listener`.
`node tool/audit-summary.mjs <report.json> [memory.jsonl …]` reads a run's report: which audits ran,
the findings by rule, by feature (leaked scopes, timers and listeners apart) and by scenario, the most
frequent ones, and with `BDD_MEMORY_LOG` files each worker's heap and renderer from its first feature
to its last.

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

**Kinds** cover the whole u2 library and EMS: one noun for every `data-u2` kind they stamp (inputs,
buttons, forms, lists, trees, tables, dialogs, tabs, sections, accordions, cards, wizards, menus,
notifications, tooltips, the filter builder, the domain controls …), generated from the libraries'
manifests into `bindings/common/u2-kinds.g.ts`; `kinds.ts` adds only what a manifest cannot know —
the Dart shell's selectors and conventions (`icon-scatter-plot`, `viewer-Grid`,
`div-section--Viewers`, `input-host-Caption`), legacy aliases, parts and gestures — and the shell's
`viewer`, `view`, `element`. A `readout` is a "name = value" status line (`[role="status"]` with a
`<label>` and an `<output>`): `value of rows readout` is the value alone. A row goes by its primary
text. Every suffix split is kept: `sequence column input` is tried as `column input` qualified
"sequence" and as `input` qualified "sequence column". `grok-bdd lint` prints how each phrase resolved.

**What the user sees comes first.** A qualifier narrows a kind's candidates by what a person sees —
the kind's own order among its label (an input's label part, a dialog's or a section's title), the
accessible name the browser associates (`<label for>`, `aria-labelledby`, `aria-label`, right after
the label), the title, the text, the placeholder and an item's key — and by `data-u2-name` only
after all of them; a Dart `name=` keeps the place its kind gives it. The first that finds something
wins: `"Search users…" typeahead` is the typeahead showing that placeholder, and a control named in
code is reached by the label a person reads. `BDD_RESOLVE_TRACE=<file>` writes a
JSON line for every phrase this order resolves by another strategy than the kind declares, and
whether both reach the same elements.

**Items of a u2 collection** (an item or option, a tree node, a table row, a tab, a wizard step, a
menu item, a breadcrumb, a chip) match by what they show, then by their key (`data-u2-item`: the
value, the tab id, the menu path `File > Open`). An item that is not rendered — a virtual list's
row scrolled away, a closed popup's option — is brought into view by u2 (`DG.U2.revealItem`) when
the phrase names its collection (`<item> in <collection>`) and nothing visible matches.

**Several matches are an error for a gesture.** A click, a hover or a check whose phrase matches
more than one visible element fails naming each candidate the way a person tells them apart — its
label, the section, dialog, form, tab or card around it, and `[name=…]` where the label differs —
instead of acting on the first (typing still goes to the first match's editor).

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
a generic kind of its own with `kind()` for a repeated structure that carries no `data-u2` (Chem's
`calculator`); markup a library kind already reads needs none — U2Demo's status lines are
`readout`s. Selectors inside `labelSelector`, `parts` and `has`-style filters are
evaluated *from the element*: `span:first-child`, not `.x > span:first-child`.

## Steps

`grok-bdd list-steps` prints every phrase the package can use and where it comes from — that
list is the reference; this is the map:

- **Gestures and outcomes on any element** (`bindings/common/steps.ts`): clicks, hovers, typing
  (`types` / `enters` = types and commits), keys, `selects`, checks, expands, drags, `fills in:`;
  `should be/become {state}`, text, value and item counts, a visible count remembered and then
  claimed `fewer`/`more … than remembered` (a search that narrows a stand-sized list), a value
  remembered and claimed again after a reload (`should (not) have the remembered value`), one choice
  a dropdown does or does not offer, and the downloads (none at all, a file by a name pattern, the
  occurrences of a text in the last one). States: visible, hidden, present,
  absent, enabled, disabled, checked, unchecked, partially checked, selected, empty, expanded,
  collapsed, focused, invalid, valid, ready, read-only, editable, busy, required — each read from
  the u2 status or the ARIA state the element uses (below).
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
  `remove button` and `checkbox` parts, typed into through `membership search`. Signing in as
  another account on the feature's page (the sharing user or a named fixture account) and back,
  the running account back first thing at feature end; the second account's notifications.
- **The workspace** (`platform/workspace.ts`): what is open (no table left, the table views, a
  table view opened with its rows, a project), a project's direct link, a file put into the user's
  files or a space for the feature (deleted at feature end), a file dropped from disk, the file
  chooser a menu opens.
- **The current table through the JS API** (`platform/data.ts`, `columns.ts`): selection and
  filter set and checked row by row (and how many selected rows the filter lets through), cells
  (every value, some value, distinct lengths, two columns equal row by row), calculated and renamed
  columns, colour coding, other open tables, links between tables, the filter panel's cards through
  its own API.
- **The top menu and its commands** (`platform/commands.ts`): a path picked by real pointer moves,
  the function call it starts awaited, the columns it added read back.
- **A package function the UI offers no entry to**, called for what it shows (`platform/functions.ts`;
  what it returns is a package test), **custom platform events**
  and the task bar's progress entries (`platform/events.ts`), the clipboard and a file chooser
  (`common/steps.ts`).
- **u2 controls through their status** (`common/u2-steps.ts`) and **geometry** (`common/geometry.ts`),
  tabled below.
- **The `viewers` tier** (below).

**The status first.** On the root of a u2 control, a state or a value claim reads what the control
reports (`DG.U2.state`) and falls back to the DOM for a key it does not report; a part, a wrapper
and a Dart widget are read from the DOM. When both answer and disagree, the claim follows the status
and prints `bdd: u2 contract: …` once per phrase and page — a finding for u2. A value failure names
what was read and the control: `the value of Stage input is "Discovery", expected "Phase II" (u2:
choice-input "Stage", enabled, valid)`.

**The visible match first.** A value or a state reads the first visible match, else the first (a tab
panel not shown keeps its own copy first in the DOM); `contain`/`have text` over several matches
reads the visible ones when there are any, else all — any of them for a positive claim, none for a
negative. A hidden match never satisfies a claim while a visible one exists. Several visible matches
keep each step's rule: the u2 steps fail naming them, a value reads the first, a flag state
(`enabled`, `invalid` …) keeps its every/any reading.

The u2 steps answer from the status alone and fail on an element that is not a u2 control root
(the error steps fall back to the error part):

| Phrase | Claims |
|---|---|
| `{element} should have the error {string}`, `… should have no error` | the validation message, exactly |
| `{element} should have {int} item(s) in total` | every item, the rows a virtual list has not rendered included (`should have {int} item(s)` counts the rendered ones) |
| `{element} should (not) list the item {string}` | an item, rendered or not, by key, else label |
| `item {string} of {element} should be {state}` | one item: selected, (un/partially) checked, expanded, collapsed, enabled, disabled, valid, invalid |
| `{element} should have the selection {string}`, `… should have nothing selected` | the selected items as a set, a comma-separated list of keys or labels |
| `{element} should be on the {string} step`, `{element} should show the {string} tab` | the current step or tab, by key, else title |
| `{element} should report the value {string}` | the raw value where the shown text formats it: as JSON when the argument parses, else as text |

The geometry claims poll the box of the one visible match of each phrase, allow 1 px and fail with
both boxes in one sentence:

| Phrase | Claims |
|---|---|
| `{element} should be left of / right of / above / below {element}` | the facing edges do not cross |
| `{element} should be inside {element}` | every edge within the other's box |
| `{element} should be inside the window` | every edge within the browser window — a dialog or popup that no ancestor clips |
| `the window is {int} by {int} px` | resizes the browser window until the test ends (the shell reset puts back the run's size) |
| `{element} should not overlap {element}` | the boxes share at most 1 px |
| `{element} should be left-aligned / top-aligned / vertically centered with {element}` | the left edges, top edges or middles agree |
| `{element} should be as tall as {element}`, `{element} should be {int} px tall` | the height — a u2 editor beside the Dart one it matches |
| `{element} should be as wide as {element}`, `{element} should be {int} px wide`, `{element} should be at least {int} px wide`, `{element} should be at most {int} px wide` | the width — an editor as wide as its Dart twin, no shorter than a floor, no wider than a cap |
| `{element} should not be clipped` | inside every ancestor that clips it, its text not cut or ellipsized |
| `{element} should have no layout problems` | the layout audit ("Audits", above) over its u2 controls finds nothing |

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
  value range and the colour scale against before; an area `painted in no color` (greys only).
- **The legend** (read from its `data-legend-*` attributes: its mode — docked, in a corner,
  collapsed to the mini icon, shown in the mini icon's tooltip, placed nowhere — its slot, its items
  and their colors, its size against before after a splitter drag, whether its box stays clear of a
  hit area, whether its items are drawn as structures or as text), the row tooltip, viewer events,
  layouts saved and loaded, sizes held and restored, and the floors: `no errors should have been
  logged`, `no error or warning balloon should have been shown`; a balloon that should have been
  shown by kind and text (`an error or warning balloon matching "<regex>"` for either kind).
- **`widgets.ts`** holds the steps first written for one viewer that a second wanted: the viewer's
  own menu, the description's place, empty plot space, range sliders, on-viewer column selectors,
  inner viewers, card readings, lassos, cross-widget drags, tabbed panels, grid pins, a right-click
  on the empty space past a grid row's last column, the columns a column picker offers (read from
  its own table and dismissed), and a Select columns list exactly or without a column.
- **`formula-lines.ts`**: the Formula Lines dialog (a line added, its range, an edit), the lines a
  viewer draws (the `formula line <title>` / `formula band <title>` areas it reports only for what
  it drew), their ranges, and the formulas of the viewer and the table naming only existing columns.
- **`filter-panel.ts`**: the cards of the filter panel, what every view's panel filters by across a
  project round trip, and a wait for the Chem cards (searching, drawing, scaffold hits).

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

**The `molecules` tier** is for a package whose features draw or read molecules wherever a host opens a sketcher (a
cell editor, a molecule input's dialog, an app's inline sketcher): `crux sketcher widget` (Crux Sketch, Chem's own
sketcher, its root named `cruxSketch`: its atoms and bonds as `atom N` and `bond N` areas; its parts as areas and as
`<part> of crux sketcher widget`, `canvas`, `actions`, `tools`, `elements`, `templates` and `label editor`; each control
a toolbar shows as a `tool <name>` area, by Crux's name for it (`tool ring.benzene`); its `ready`, `pending`, `smiles`,
`atoms`, `bonds`, `mode`, `selected atoms`, `selected bonds`, `tool`, `query`, `empty` and `changes` readings; its
gestures settle on its `isRenderPending` and `onRendered`), Crux's controls by their test ids (`crux canvas`, `crux single bond tool`, `crux benzene tool`,
`crux nitrogen tool`, `crux R-group tool`, `crux R1 button`, `crux clear button`, …), `the "<reading>" reading of
<widget> should be the molecule "<SMILES>"` and `… the molecule in row N of "<column>" column`, `the molecule in row N
of "<column>" column should be "<SMILES>"` (each read by RDKit, any notation), and `the semantic types of the current
table are detected`. It needs Chem on the stand, so a package other than Chem gates on it (`Given the "Chem" package
is installed`); a feature that draws on Crux pins it (`the molecule sketcher is "Crux"`) and is `@sketcher-controls`.

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

A widget that is no viewer and lives outside one (a molecule sketcher in a dialog) takes part the
same way: its own status provider, registered on itself, makes it a widget the platform knows, and
the viewers tier's area and reading steps find it as the nearest such widget up from the element a
phrase names. Chem's Crux sketcher reports its atoms and bonds as `atom N` and `bond N` areas and
readings such as `smiles`, `atoms` and `pending`, and a `get isRenderPending()` its gestures settle on:

```gherkin
When user clicks on crux single bond tool
And user drags the "atom 1" area of crux sketcher widget by 60 pixels to the right
Then the "atoms" reading of crux sketcher widget should be 3
```

## Guides: a scenario as a how-to video

A scenario that answers "how do I …" is also its own demonstration. `grok-bdd guide
features/guides/<name>.feature` compiles it, runs it on one worker in guide mode and renders each
scenario into `guides/<feature slug>/<scenario slug>/`:

- `guide.mp4` — the pointer travels to every element a step acts on, the element is lit (the rest
  of the page dimmed) and rests under the pointer before the click lands; every press is marked
  where the page received it, under the pointer's tip — a yellow dot and ring for the left button,
  a green one for the right, twice for a double-click; the pointer never skips, each movement
  starting where the previous action ended (the renderer refuses to make a video in which it
  skips); an icon-sized target
  (28 px or less each way: an icon, a checkbox) is zoomed into first, anything larger is clicked
  where it is; a gesture on a hit area of a viewer or a widget (a grid cell, an atom of a sketcher)
  lights that area, not the viewer round it; the page after the step is revealed, and a caption above the page (clear of a
  player's timeline) reads the step as an instruction ("Click on Open local file icon in browse
  toolbar"), or says what a `# caption: <text>` comment on the lines right above the step says — for a
  step whose text is for the test (an atom's index, a drag in pixels, a reading), so the video speaks as
  a person would ("Click the acid's carbon", "The centre flips from R to S"); a check with a caption is
  shown even where `HIDDEN_CHECKS` would leave it out; a caption too long for the strip is set smaller, then on two lines, and the lines of a
  pasted text (`\n` in the step) read as values separated by commas. A choice is shown being made: a native `<select>` opens its list, the option is typed
  so the list highlights it, Enter takes it; a column selector opens its picker, the name is typed
  into the search short of its last letter (a complete unique name is taken on the spot) and the
  row it leaves is clicked. A `Then` step shows what it checked with a check mark —
  when a person could see it (a dialog, a column, a row count, a value, a legend item's color);
  the checks a test needs and a person does not (error and balloon floors, server state, viewer
  readings and pixels, "than before" claims, property bags, widget counts, task-bar and command
  bookkeeping, values matched against a regular expression) are left out, by the `HIDDEN_CHECKS`
  patterns in `src/runtime/guide.ts` — so a guide names a new column exactly, not by a pattern;
- `step-NN.png` — the lit picture of every step, and `steps.md` — the numbered steps with those
  pictures, ready to paste into a reply (a filming with fewer steps removes the pictures it no
  longer has);
- `audit.png` and `audit.json` — every press as the video shows it, beside the same picture with
  ticks aimed at where it landed, and how far the mark's centre and the pointer's tip are from it
  (the renderer prints a press more than 2 px off, or outside the element its stop lit);
- with `--gif` also `guide.gif` and `guide-thumb.png`, the docs' own pair (the GIF's palette
  always holds the press colours).

A set-up step can frame the video to a part of the page — `guideFrame(page, box, margin)` from
`@datagrok-libraries/bdd/runtime` (Chem's `the Crux sketcher is open on …` frames its sketcher's
dialog, so a guide about the sketcher shows the sketcher alone): the renderer cuts every picture to
that box, and moves every point and box into it.

Guide mode (`BDD_GUIDE=<dir>`, set by the command) records at the step: the page before and after
it (`BDD_GUIDE_SETTLE`, 500 ms by default, lets a dialog or a balloon finish appearing), the last
element the step located, and what the pointer did — as the page saw it: capture listeners log
every real press, release and move, whatever sent it (a locator's own click never goes through
`page.mouse`), each press with its button, its place and the element under it, and a run of moves
with no button held as the one point it came to rest. A step's presses and moves go to the stop
they were made at (below). A move with a button held is a
drag: the page is pictured along the way (`NN-dragK.png`, at most eight per step), so the video
shows what the drag draws — a selection box, an annotation region — growing under the pointer,
and the step's still shows it complete at the release point. Tests know nothing of it: without
the variable no line of it runs. The viewport is 1080p (1920×1080, so the top menu keeps every group on the bar; `BDD_GUIDE_VIEWPORT=<w>x<h>`
for another) so the video reads without zooming every step — the video itself, caption strip included, is 1600×1000, the frames downsampled from the page (`VIDEO_W`, `VIDEO_H` in the renderer) — and the shell is the full one (simple
mode off), as a person has it — filmed or in a plain run: every `@guide` scenario carries `And
simple mode is off` right after the login (the compiler refuses one without it), and the step
puts simple mode back at feature end. Every step is in the video except the login and that shell
step (`guide.silent`), other set-up a person does not take (a pinned setting such as `the molecule
sketcher is …`, a package's own readiness wait, which calls `silent` from
`@datagrok-libraries/bdd/runtime`), a wait the `HIDDEN_CHECKS` list names (`the package autostarts
have completed`, `… should have finished updating`), whether it is written as a `Given` or a
`Then`, and a step that neither
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
Guides live under `features/guides/` and run with the rest of the suite, in the full shell: an
answer that stops being true fails a test.

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

`npm run build` compiles `src/`, `bindings/` and the Playwright config to `dist/`, which is what a
package's `grok-bdd` loads: while any of those sources is newer than its build (a pull, an edit),
every command but `init` and `link` stops and names the file, rather than failing a sound feature on
a kind or a step only the sources have. `npm run
test:unit` runs the engine tests (nouns, compile, project, init, failure, vocabulary, audit …) and
the ones that drive a static page in the library's Chromium (locate, addressing, u2-state, layout …;
each skips itself where none is installed); the library is a project itself (`features/platform`):
`npm test` builds, checks the vocabulary, drift-checks the project and runs the unit tests — what the libraries CI runs on Node 18,
which is why `@playwright/test` is pinned to the last minor that runs there — and `npm run
test:suite` runs it against a stand. Translating a hand-written spec into a feature, and proving
the feature tests what it claims, is the `/bdd-translate` skill.

The u2 and EMS nouns are generated: after any change to `libraries/u2/manifest.json` or
`libraries/ems/manifest.json`, `npm run vocab` rewrites `bindings/common/u2-kinds.g.ts`; `npm run
vocab:check`, `npm test` and the vocabulary unit test fail while it is older than a manifest. What a
manifest cannot say (a Dart selector, a legacy alias, a strategy order a suite relies on) goes into
the `OVERRIDES` of `bindings/common/kinds.ts`; a compat test holds every key, selector, strategy,
part and Dart name of the hand-written table the generated one replaced
(`tests/fixtures/kinds-before.json`), but for the narrowings it lists (`REPLACED`).

# @datagrok-libraries/bdd — behavioral automation (Gherkin → Playwright)

Feature files bound to the u2/platform vocabulary, compiled deterministically into Playwright specs
that packages own. Design record: `core/docs/features/ui2/automation/BRAINSTORM.md` (rulings by the
lead: standard Gherkin; committed, drift-gated codegen; u2-centred; overridable composition;
generic kinds; reserved platform names; packages own their tests; one page per feature; a global
CLI; NOT wired into `grok test`; §8 there is the current state, gaps and next steps). The contract
it consumes and the u2 findings: `core/docs/features/ui2/AUTOMATION.md`; the lane map:
`core/docs/features/ui2/TESTING.md`.
Authoring guide: `README.md` here.

## Layout

```
bin/grok-bdd.js    launcher: prefers the package-local install (its specs import that copy), else this one; runs dist/src/cli.js
src/project.ts     a project = a dir with features/ (a package's bdd/, or the library root); config bdd.config.json {tiers}; binding sources
src/discover.ts    imports binding modules (library: this build's bindings/; project: .ts through tsx) and maps step fn → export name
src/gherkin.ts     @cucumber/gherkin → flat model (And/But resolved, outlines expanded, rule backgrounds folded)
src/registry.ts    Given/When/Then/Step, element(), context(), alias(), kind(), dataset(), defineParameterType(); one global registry
src/nouns.ts       phrase (+ context) → NounRef (pure; shared by compiler and runtime); parts of elements AND kinds
src/match.ts       cucumber-expressions matching + specificity
src/compile.ts     FeatureModel → *.test.ts (feature(test) session, test() per scenario, test.step per Gherkin step, enter/leave)
src/init.ts        `grok-bdd init`: scaffold bdd/ in a package (never overwrites; merges .vscode settings, .gitignore, package.json)
src/cli.ts         init | compile [--check] | lint | list-steps | run [playwright args]  (run = check + Playwright with playwright.config)
src/runtime/       args (el/ds/enter/leave), locate, gestures, assertions, harness (feature session + resetShell), global-setup, index
src/index.ts       what package bindings import from '@datagrok-libraries/bdd'
bindings/common/   parameter-types, kinds (ALL u2 data-u2 kinds + Dart conventions), steps, session — base, always loaded
bindings/platform/ the shell: elements (reserved names), datasets, steps — base, always loaded
bindings/tiers/<t>/ opt-in tiers (`viewers`); a project names them in bdd.config.json
features/ generated/ bdd.config.json   the library's own project (platform smoke feature, tier viewers)
playwright.config.ts  the one config every project runs with (BDD_ROOT → testDir/outputDir/storageState under the project)
dist/              `npm run build` output (tsc, ESM, .d.ts, source maps) — what `exports` and the launcher use; gitignored
tests/             node:test via tsx (nouns, compile, project, init)
```

A package project: `<pkg>/bdd/{package.json {"type":"module"}, bdd.config.json, features/, bindings/,
generated/}`; sample `packages/U2Demo/bdd`.

## Rules that must not regress

- **Everything a spec imports from the library goes through `dist/` by package subpath**
  (`@datagrok-libraries/bdd/runtime`, `…/bindings/common/steps`); project bindings are imported
  relative (`../bindings/steps.js` → Playwright maps to `.ts`). One registry instance per process is
  what makes context/element registration visible to the runtime — never mix `src/` and `dist/`
  imports in one run. The launcher resolves the package-local copy first for the same reason.
- **ESM everywhere.** The library is `"type": "module"`; a package's `bdd/package.json` marks its
  specs and bindings as ESM too — a CJS spec would `require()` the ESM dist and fail on Node 20.
  Project `.ts` bindings load at compile time through `tsx/esm/api` `register()`; at run time
  Playwright transpiles them.
- **One `@playwright/test` per run.** Specs import it from the package, the library's runtime from
  its own real path; on a junction-linked checkout the package gets a junction to the library's copy
  (`node_modules/@playwright/test`), on an npm install the peer dependency is shared. The harness
  receives `test` from the spec and never imports it.
- **One page per feature** (`src/runtime/harness.ts`): `feature(test)` registers `afterEach`
  (leave the context, `resetShell`) and `afterAll` (close the context); the page is created inside
  the first test (`session.page(browser)`), so Playwright merges the project's context options
  (storage state, viewport) and records traces/screenshots for it. The generated spec calls
  `test()` itself so reports point at the spec line, not the harness.
- **Names are global and unique; platform names are reserved.** `element()` throws on a name or
  alias already registered; a context's `element()` throws on a global name too. App vocabulary goes
  on a `context(name, def)`.
- **Context switching is compile-time tracked and runtime explicit.** Step meta `enters: '<ctx>'`
  makes the compiler validate later phrases against that context and emit `enter(page, '<ctx>')`;
  the runtime keeps the current context per page; `resetShell` leaves it. Every runtime parse goes
  through `refOf(page, target)` (locate.ts) — a parse without the page loses the context.
- **Noun resolution order** (`src/nouns.ts`): whole phrase registered (context names first, then
  global) → ordinal → split at the first scope word outside quotes (inner within outer, recursive;
  `of` names a part of a registered element or of a kind) → registered element, else generic kind by
  suffix with EVERY matching suffix kept, longest first.
- **Context-first lookup**: a context-local name or a generic kind without an explicit scope is
  searched inside the context root first, then on the whole page (portaled dialogs, notifications).
- **Kind qualifier strategies** run in the kind's `match` order, first hit wins, union (`.or`) when
  nothing matches (negative assertions still get a locator). Kinds cover every u2 `data-u2` value
  (`grep -rhoE "dataset\.u2 = '[^']+'" libraries/u2/src`) AND Dart `name=` conventions. Plain u2
  `button()`, toolbar buttons and tab headers carry NO `data-u2` — those kinds match by tag/role +
  text. Qualifiers are lowercased; text matching is case-insensitive (`exactText()` regex).
- **Owner edge**: nothing inside the outer → retry inside `[data-u2-owner="<outer's name>"]`.
- **Playwright scopes inner selectors to the element**: a `labelSelector`, a part, a `has:` filter
  is evaluated from the outer element, so `.u2demo-status > span:first-child` never matches the
  first span OF a `.u2demo-status` — write `span:first-child`. Cost a full failed run (2026-09-03).
- **`page.evaluate` must return nothing DOM/Dart-shaped**: returning `grok.shell.addView(...)`
  fails with "Cannot serialize result: object reference chain is too long" — wrap in a block and
  return undefined.
- **States** (`src/runtime/assertions.ts`): `selected` = `loc.and('[aria-selected=true], [aria-pressed=true],
  [aria-checked=true], [aria-current]')`; `expanded`/`collapsed` read the element's own
  `aria-expanded` or the first one inside (section/accordion/category headers); `visible`/`hidden`
  over several matches = any/none (`filter({visible: true})`) — a bare `toBeHidden` on a multi-match
  locator is a strict-mode error. `expectCount` prefers `tbody tr` (no header row).
- **Gestures** (`src/runtime/gestures.ts`): `editorOf` = the element itself when editable, else
  the first of `EDITOR` (inputs, selects, textareas, contenteditables), else the
  `[data-u2-part="editor"]` trigger (icon/function/columns pickers), else the element (a slider
  handle, a list). `select` = native `<select>` first; else click the editor, ArrowDown when it is
  a `role=combobox` (u2 comboboxes/typeaheads open on a keystroke, not a click), then the option
  by whole text → primary-text part (`OPTION_LABEL`) → `title`/`aria-label` (icon cells) → substring.
  `setExpanded` clicks the twistie of a tree row (a row click only selects), else the
  `aria-expanded` control. `fillIn` recognises switches through `role=switch` inside the editor.
- **Step specificity**: fewer parameters, then more literal characters; ties are ambiguous errors.
  `{key}` has no spaces so `user presses {key} in {element}` cannot lose to `user presses {key}`.
- **Codegen** emits names, never selectors; `\n` line endings, no timestamps, EOL-normalized drift
  check; orphaned generated files are removed on compile and reported on `--check`.
- **Session readiness**: `user is logged in` skips navigation when the page is already in the shell
  (chained scenarios), calls `grok.shell.closeAll()` and waits for the Home view to be current
  (`grok.shell.v.type === 'datagrok'`) — closeAll re-adds Home asynchronously.
- Do not reinvent what `@datagrok-libraries/test/src/playwright/*` has (storage-state login,
  `openTableFromFile`, base config): import it (compiled `.js` + `.d.ts`, `.js` suffix).

## Environment facts that cost time to learn

- Dev setup on this checkout: the library is junction-linked into U2Demo
  (`packages/U2Demo/node_modules/@datagrok-libraries/bdd` → `libraries/bdd`, and
  `node_modules/@playwright/test` → the library's copy); create junctions with PowerShell
  `New-Item -ItemType Junction` — Git Bash mangles `mklink /J`. Run `npm run build` in the library
  after engine or binding changes; the U2Demo project is driven with
  `node ../../libraries/bdd/bin/grok-bdd.js <command>` from `packages/U2Demo` (a global install /
  `npm link` gives plain `grok-bdd`).
- Stand: `http://localhost:8888` (nginx over datlas `:8082`), admin/admin, developer key `admin`.
  When every login "hangs", the hand-started Postgres container (`affectionate_einstein`, no restart
  policy) died with a Docker restart: `docker start affectionate_einstein`.
- To publish U2Demo: `cd packages/U2Demo && npx webpack && grok publish localhost --key admin
  --skip-check`. Webpack and `grok check` ignore `bdd/`; U2Demo un-ignores `src/raw.d.ts`.
- U2Demo/u2 builds: `@datagrok-libraries/u2` is NOT on npm; `datagrok-api` in `libraries/u2` and
  `packages/U2Demo` node_modules must be declaration-only copies of `public/js-api`; `grok api`
  rewrites `package-api.ts`/`package.g.ts` with LF endings (content-identical; `git checkout --`).
- In the U2Demo app view the platform `.d4-toolbox` is present and visible but empty; the app's
  navigation is a u2 splitter panel of its own. U2Demo's `bdd/bindings/demo.ts` makes the
  sub-demo content pane (`.u2demo-content`) the `U2 Demo` context so page trees/lists win over
  the navigation tree (`demo navigation` by name), registers a package kind `readout` for the
  pages' "name = value" lines, and opens a sub-demo with `user opens the "<label>" demo page`
  (cold page: `page.goto`; warm page: `grok.functions.call('U2Demo:u2DemoApp', {path})` +
  `grok.shell.addView`, so tables opened by earlier steps survive — a `goto` reloads the client;
  the shared helper is `openSubDemo(page, route, ready)`, which the MSA workbench step uses too).
  One feature per sub-demo lives under `bdd/features/demo/<group>/`.
- **Live DOM outline probe**: before writing phrases for a page, dump what the engine will see.
  Recipe (scratch dir): `createRequire('<libraries/bdd>/package.json')('@playwright/test').chromium`,
  `newContext({storageState: '<pkg>/bdd/e2e/.auth.json'})`, goto the page, then `page.evaluate` a
  walker over the root printing tag, classes, `data-u2*`, `role`, `aria-*`, `name`/`placeholder`,
  input values and own text per element (mark `offsetParent === null` as hidden); for popups dump
  `[data-u2="menu"], [data-u2="dialog"], [data-u2="tooltip"], [data-u2="notify"], [data-u2="tour"],
  .d4-dialog, .d4-balloon` after the opening action. Pass the walker as a string expression that
  CALLS the function (`page.evaluate(\`(${fn})(${JSON.stringify(args)})\`)`) — a bare function
  string is evaluated, not invoked. u2 popups are portaled under `body > .u2-overlay` with
  `data-u2-owner` = the nearest NAMED ancestor of the trigger (the app root, "U2 Demo", when the
  input has no name), so the owner edge rarely helps inside the demo — the whole-page fallback does.
- Facts the probe established (2026-09-03): u2 ChoiceInput is a native `<select>` with a leading
  empty option; BoolInput `switch` is `span[role=switch][aria-checked]` (no input); NumberInput's
  editor part wraps `span > input[type=text]`; TextInput `search` wraps the input in the editor
  part; icon/function/columns pickers are `[data-u2-part=editor][role=button][aria-haspopup]`;
  Combobox has no input root (`[data-u2=combobox] > input[role=combobox]`) and opens on input or
  ArrowDown; MultiSelect's field is `[role=combobox]` opening a `.u2-multi-select-popup` with
  `[role=option]` rows and a `Select all` `[role=button]`; ButtonGroup items are `button`s
  (`role=radio` + `aria-checked` in single-toggle mode, `aria-pressed` in multi); toolbar toggles
  carry `aria-pressed`; menu items are `.u2-menu-item[role=menuitem]` with `.u2-menu-label` and
  `.u2-menu-shortcut`, `aria-disabled`, submenu `aria-haspopup`; menu-bar items are
  `button.u2-menu-bar-item[role=menuitem]`; breadcrumbs are `button.u2-breadcrumbs-item` +
  `span.u2-breadcrumbs-current[aria-current=page]`; accordion panes `.u2-accordion-pane` with a
  `[role=button][aria-expanded]` header and `.u2-accordion-title`; Section header
  `.u2-section-header[role=button][aria-expanded]`; Wizard steps `li.u2-wizard-step` with
  `.u2-wizard-title` and `aria-current=step`, buttons BACK/NEXT/FINISH; Tour overlay
  `[data-u2=tour]` + sibling `.u2-tour-popup[role=dialog]` with SKIP/NEXT/DONE and a "1 / 4"
  counter, `onFinish` gets `done`/`skipped`; notify balloons `.u2-notify[role=status|alert]` with
  aria-label Close/Copy icon buttons; BasicTable is a real `<table>` (`tbody tr[aria-selected]`);
  VirtualGrid/VirtualList are `[role=listbox]` with `[role=option]` cells/rows (`.u2-list-row`,
  `.u2-grid-cell[title]`); VirtualTree rows `[role=treeitem][aria-expanded]` with
  `.u2-tree-twistie`/`.u2-tree-label` (the twistie toggles, the row selects); PropertyGrid rows
  `.u2-propgrid-row` (`.u2-propgrid-name` + inline inputs) under `.u2-propgrid-category[role=button]
  [aria-expanded]`; Card `[data-u2=card]` with `.u2-card-title`, `role=button` when clickable,
  `aria-pressed` when selectable; StatCard `.u2-stat-label`/`.u2-stat-value`; ProgressBar
  `[role=progressbar]` + `.u2-progress-percent` + `.u2-progress-description`; MessageInput editor
  `.u2-msg-editor[contenteditable][role=textbox]` + `button.u2-msg-send`; the Dart `ui.dialog`
  is `.d4-dialog[role=dialog][name="dialog-<title>"]` with `.d4-dialog-title` and
  `button.ui-btn[name="button-OK"]`; a bridged u2 input inside it keeps `data-u2="text-input"` on
  the `.ui-input-root` with a `.ui-input-label`; the shell status bar is `.layout-status-bar` with
  the view's panels in `.d4-view-status-panel`; a platform "Debugging packages" `.d4-balloon.warning`
  is up after every reload (text checks on `notification` are any-of, so it does not interfere).
- Node prints "Importing JSON modules is experimental" whenever tsx registers; the launcher
  filters that one warning.
- **VS Code Cucumber extension** (`cucumberopen.cucumber-official` 1.11.0, language-server 1.7.0;
  the lead also has `alexkrechik.cucumberautocomplete`, unconfigured). It parses
  `defineParameterType({name, regexp: /…/})` calls from the glue, so `cucumber.parameterTypes` in
  settings holds only `state` (its regexp is built from a list in code); listing the others again
  yields "already a parameter type" errors. Settings files live at the core root (gitignored), in
  `packages/U2Demo/.vscode` and here. **A settings JSON written through a Bash heredoc lost the
  `\\` of `\\w`** (2026-09-03): VS Code's lenient parser kept `[w -]+?`, so `{viewer}` and
  `{dataset}` steps showed as undefined ("scatter plot viewer should be added …") while the compiler
  resolved them — write JSON with the Write tool and validate with `json.load`. To see exactly what
  the extension resolves: `npm i @cucumber/language-service` in a scratch dir, `new
  WasmParserAdapter('<pkg>/dist')` (the Node adapter has no TypeScript grammar), sources
  `{languageName: 'tsx', uri, content}`, `new ExpressionBuilder(adapter).build(sources,
  parameterTypesFromSettings)`, then `expressionLinks[].expression.match(stepText)`.
- Bash tool here: cwd persists across calls and very long heredocs fail to parse — use absolute
  paths and the Write tool for whole files.
- **u2 finding, not a harness bug (2026-09-03, unreported)**: in a `funcForm` a number field
  typed key by key starting with "-" ends up EMPTY — the transient "-" parses to null, the form
  writes null into the FuncCall, and the param's `onChanged` echo writes null back into the input
  after the next keystroke landed (probe: `scratchpad/probe/numinput.mjs`; "7" sticks, "-7" does
  not, `fill('-5')` does). A human loses the sign the same way. The U2Demo features use positive
  numbers meanwhile; do not add a `fill` fallback to `typeInto` to hide it.
- Two phrase traps met while writing the U2Demo features: a qualifier starting with an ordinal
  word (`First name input` → "first" + `name input`) needs quotes (`"First name" input`), and a
  label equal to a kind name (`Columns input`) now resolves to the labelled input first (nouns.ts
  puts the whole-kind reading last).

## Conventions

TS strict, 2-space, single quotes, `.js` suffixes on relative imports (NodeNext), comments near
zero. Nothing is committed without the lead's order; stage by explicit path — the tree carries
other sessions' work.

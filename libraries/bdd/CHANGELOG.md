# BDD library changelog

## v.next

* Second account: Fixed a scenario that stopped while signed in as another account running every scenario after it as that account: a journey scenario that stopped, and a test that failed or is a known failure, sign the running account back in before the next one starts (`atScenarioStop`)
* Locator: Fixed a gesture on a gallery card or link landing before the gallery had redrawn its search or refresh: a phrase naming a gallery or a card waits for the gallery grid's `aria-busy` to clear
* Vocabulary: Added `user double-clicks on {element} and the project opens` (`grok.events.onProjectOpened` on a 120 s budget, kept in a `@known-failure` scenario; an error balloon or the Data loading error dialog ends it as a failure) and `user double-clicks on {element} and the query it runs completes`
* Server fixtures: A killed run's folder of a file put into the user's files (`BDDx{time}/data.csv`) is swept by its family, dated by its `{time}`
* Harness: Fixed a gallery view mode reaching the next feature on the page: the shell reset clears the galleries' view modes with the panes' expanded state (a feature that switched the Scripts view to cards left the next one without its script links)
* Vocabulary: Added `the project named {string} should be tagged {string} on the server`, `the query {string} on the server should (not) read {string}` and `the sharing user has no grant on the {string} connection, now and when the feature ends`
* Cleanup: Fixed the picture of a project deleted through the UI staying on the server: the harness records the pictures a page uploads during a feature and deletes them at its end
* Harness: Fixed a failed scenario leaving its dialog over the next one when Escape does not close it (an editor that keeps Escape, such as the formula editor of Add New Column): the shell reset and a failed journey scenario click the topmost dialog's close icon, so a failure no longer leaves its dialog over the next scenario
* Second account: Fixed a page reloaded after `the sharing user has no notifications` coming up signed in as the second account: minting its session took the login's Set-Cookie, which the page then sent on its next load
* Context menus: Fixed a right-click on a tree row under a sticky group header (a root space under Spaces once the Browse tree is scrolled a row) landing on the header three times: the row is scrolled to the middle of its panel; a failed menu pick names what the right-click reached
* Harness: A failed journey scenario closes only the dialogs it opened, not one the Background opened under a menu Escape leaves open
* Balloons: Fixed a failed balloon claim always saying "none": it lists the balloons the poll saw
* Locator: Fixed an input phrase matching a host input whose nested input carries the label (a multi-choice by its items' labels), so a claim no longer passes on the host; a u2 bridge still goes by its inner Dart root's label
* Kinds: Extended the `text input` kind to a map input's named row editor ("Pattern 1" of a conditional coloring)
* Kinds: Added a `color bar` part to the `item` kind, the bar a category of the categorical color coding editor opens its picker from
* Vocabulary: Added `the server has logged activity on the {string} script` (the Activity pane counts once, when it renders: a fresh view opened before the server logged a run showed 0 for good)
* Vocabulary: Added `the text of {string} column should not be color-coded` and `user scrolls {element} into view` (the wheel over its panel until the whole element shows, for one at the start of a list)
* Server fixtures: Fixed the stale-fixture sweep missing a `{time}` with no dash before it, as a script a formula calls by name has (`BddDouble{time}`)
* Layout audit: `label-column` measures a cell of a form's columns (u2 `Form.addColumns`, `data-u2-part="columns"`) from its own column, so two aligned columns are no longer reported as one misaligned one
* Vocabulary: Added the u2 `status-line` kind ("status line")
* Vocabulary: `the top menu should list {string} as disabled` / `as enabled` — the platform greys out a command whose function needs a column the table lacks (#4226)
* Kinds: A plain ARIA list (`role="list"`) is a list and its `role="listitem"` rows are items, so "first item in annotations list" reaches them

* Geometry: Added width claims — `{element} should be as wide as {element}`, `… should be {int} px wide`, `… should be at least {int} px wide` and `… should be at most {int} px wide` — so an editor too short next to its Dart twin fails a feature
* Geometry: Added `{element} should be inside the window` (a dialog or popup no ancestor clips) and `the window is {int} by {int} px`, which resizes the window until the test ends
* Layout audit: Inputs inside another input's editor (a columns map's key rows) no longer count as rows of the form in the label-column rule
* Runner: Added `BDD_VIEWPORT=<w>x<h>`, which runs a suite at another screen size than 1920×1080, so the layout audit checks it there
* Vocabulary: The u2 and EMS nouns are generated from the manifests' kind catalogs (`tool/gen-u2-kinds.mjs` → `bindings/common/u2-kinds.g.ts`, `npm run vocab`; `kindDefsFrom`/`registerU2Kinds`); `kinds.ts` keeps only what u2 cannot know — Dart selectors and names, legacy aliases and parts, strategy pins; container parts are anchored at their component's root; `npm test` fails on a stale vocabulary, and a snapshot of the kind table guards compatibility
* Kinds: Added nouns for the new u2 and EMS kinds: access grid, chips input, filter builder, filter query, filter group, filter row, read-only field, pick input, ref input, the domain and manifest kinds, layout panel, and the vertical and horizontal stacks
* Kinds: Added the generic `readout` kind: a `[role="status"]` "name = value" line (`<label>` + `<output>`), where "value of <name> readout" is the value alone
* Kinds: An entity card is found by its title
* Kinds: Fixed a u2 input that wraps a Dart input matching "<label> input" twice
* Addressing: Element phrases resolve by what the user sees first — the label, then the label the browser associates (`<label for>`, `aria-labelledby`, `aria-label`), then the title, text and placeholder; `data-u2-name` is the fallback, and Dart names keep their place
* Addressing: Items of u2 collections (item, tree node, table row, tab, wizard step, menu item, breadcrumb, chip) resolve by their `data-u2-item` key after what the user sees, and an item a gesture names is revealed through `DG.U2.revealItem` when it is not rendered; the chips of a chips input and the rows of a data table can be reached
* Addressing: An ambiguous gesture fails naming each candidate by its label and section, also under tsx
* Addressing: Added `BDD_RESOLVE_TRACE=<file>`, which logs every phrase the new order resolves differently and whether it reached the same element
* State claims read the u2 status: on a u2 control root, enabled/disabled, checked, selected, expanded, focused, invalid/valid, ready, empty and values come from `DG.U2.state`, with the DOM as the fallback per key; a status/DOM disagreement prints `bdd: u2 contract: …`
* Fixed value, state and text claims reading a hidden first match (a tab panel not shown) while a visible one exists
* Added the states read-only, editable, busy and required
* Added u2 phrases for a control's error, its total items, listed items, one item's state, the selection, the current tab or step, and the raw value; several visible matches fail naming the candidates
* A value failure names what was read and the control: `the value of Stage input is "Discovery", expected "Phase II" (u2: choice-input "Stage", enabled, valid)`
* Added geometry claims: `{element} should be left of / right of / above / below / inside {element}`, `should not overlap`, `should be left-aligned / top-aligned / vertically centered with`, `should be as tall as`, `should be {int} px tall`, `should not be clipped` and `should have no layout problems`; they allow 1 px, poll, and fail with both boxes in one sentence
* Added the layout-audit engine (`auditLayout` in `src/runtime/layout.ts`): over the visible u2 controls it finds zero-size, offscreen, overlap, clipped and truncated text, overflow, label/editor misalignment and form label-column problems; it compares only within one layer, waits on `aria-busy`, takes an allowlist and uses no screenshots
* Added automatic audits after every passing scenario on a u2 page: the layout audit and the console-error floor on the final state before the shell reset, and a leak check after it against a baseline taken when the scenario started (`DG.U2.debug`)
* Added an `audits` block to `bdd.config.json`: layout, errors and leaks set to fail, report or off, with allowlists whose entries each need a reason; `BDD_AUDIT` and `BDD_AUDIT_LAYOUT|ERRORS|LEAKS` override it; the defaults are layout and errors report, leaks off; u2 tracking (`__U2_DEBUG__`) is switched on only in pages of projects that audit leaks
* Added `tool/audit-summary.mjs`: findings by rule, by feature (leaked scopes, timers and listeners) and by scenario, the most frequent ones, and per-worker heap and renderer readings from `BDD_MEMORY_LOG`

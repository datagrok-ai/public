# BDD library changelog

## v.next

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

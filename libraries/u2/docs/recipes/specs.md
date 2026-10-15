# Recipe: dg-ui/1 specs — binding, events, plugins

A spec is JSON: `{"$schema": "dg-ui/1", "components": [...], "root": {...}}`. Each node is
`{tag, name?, props?, bind?, on?, children?}`; `manifest.json` lists every tag, its props and
events. This recipe is the composition-level rule set a spec author (or an LLM) needs on top of it.

## Binding: every prop is live

```json
{"tag": "u2-button", "name": "save", "bind": {"text": "$.saveCaption", "visible": "$.canSave"}}
```

- Every declared prop accepts a `bind` entry and **follows its source in place**: the node is
  never rebuilt, so focus, caret, scroll and unbound typed values survive. A prop marked
  `bindable` takes the signal itself (`twoWay` also writes edits back); every other prop is written
  through the component's setter on each change.
- A prop marked **`structural`** (an input's `inline`, a splitter's `direction`, a form's `layout`,
  a source's `func`) can only change by a rebuild: a bind reads it **once** and the console says so.
  Change it in the designer (a patch) or in the document, not through a bind.
- Every visual tag has a live `visible` prop. Child props a container reads off its children —
  a tab's, accordion pane's or wizard step's `title` — are live too; their `icon`, a step's `id`
  and `done` are structural.
- A bind on an undeclared prop is an error; a node may bind to a node declared later in the
  document (resolved when the render pass ends, still without a rebuild).
- A malformed node — a `bind` value that is not a `$.` path string, `children` that is not an
  array of nodes — becomes a visible placeholder naming the mistake; the rest renders.

## Events

`on` maps an event to a command: `"cmd:save"` (a context command, then a platform function),
`"cmd:Pkg:fn"`, `"cmd:grid.refresh"` (a named node's function), or
`{"cmd": "save", "args": {"id": "$.list.selected"}}`. Specs never carry code.

```json
{"tag": "u2-wizard", "name": "importer",
 "on": {"finish": "cmd:importAll", "cancel": "cmd:closeImport", "step": "cmd:trackStep"},
 "children": [
   {"tag": "u2-panel", "props": {"title": "Source", "id": "source"}, "children": [...]},
   {"tag": "u2-panel", "props": {"title": "Options", "id": "options"}, "children": [...]}]}
```

`u2-wizard` fires `finish`, `cancel` and `step`; bind a readout or a source param to
`$.importer.currentStep` (the open step's id). Gates, commits and step actions are code options
(see [forms-and-details](forms-and-details.md#multi-step-flows)).

## Companion libraries register through a plugin

A spec naming `u2-domain-*` tags needs ems installed into the registry it renders over, once,
before the first render:

```ts
import {installEms} from '@datagrok-libraries/ems/src/dg';
installEms();                      // the global registry; idempotent; answers the uninstall
```

`installEms()` is `registry.use(emsPlugin)`: a `U2Plugin` brings its tags and its spec ambient
(the ems unit of work) together, once per name. A spec naming a known plugin's tag before it is
installed renders a placeholder that says which call is missing.

# u2/dg — the platform layer

The only part of u2 that imports `datagrok-api` — eslint enforces it (D3 in the u2 plan,
`core/docs/features/ui2/PLAN.md`). Everything here bridges the platform-free core onto Datagrok
lifecycle, events, and value editors.

| Export | What it does |
|---|---|
| `host(component, closeIn?)` | Component → `DG.Widget`; disposal joins `Widget.subs` |
| `toSignal` / `toObservable` | rxjs ↔ signals, both directions |
| `leakReport()` | live u2 scopes vs registered widgets |
| `asDartInput(input, {dataType?})` | u2 `Input` → `DG.JsInputBase` (a `DG.InputBase` everywhere) |
| `columnInput(label, table, {filter?})` / `tableInput(label)` | platform-bound `ChoiceInput`s |
| `propertyForm(props, source, options?)` / `objectForm(source, options?)` | a whole form out of `Property` metadata |

Import path today: `@datagrok-libraries/u2/src/dg/index.js` (the `@datagrok-libraries/u2/dg`
subpath needs an `exports` map in `package.json`).

## Registering a u2 input as a platform value editor

A function with `role: 'valueEditor'` returning a `DG.InputBase` becomes the editor the platform
picks for every property with the matching `propertyType`/`semType` — property panels, function
dialogs, `InputForm`, `Dialog.add()`:

```ts
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {TextInput} from '@datagrok-libraries/u2';
import {asDartInput} from '@datagrok-libraries/u2/src/dg/index.js';

export class PackageFunctions {
  @grok.decorators.func({
    meta: {propertyType: 'string', semType: 'Barcode', role: 'valueEditor'},
    outputs: [{type: 'object', name: 'result'}],
  })
  static barcodeInput(): DG.InputBase {
    const input = new TextInput({label: 'Barcode', placeholder: 'SC-000000'});
    input.addValidator((v) => /^SC-\d{6}$/.test(v) ? null : 'Expected SC-NNNNNN');
    return asDartInput(input);
  }
}
```

Or ad hoc, without registration:

```ts
const barcode = new TextInput({label: 'Barcode'});
const dlg = ui.dialog('Register').add(asDartInput(barcode).root).onOK(() => save(barcode.value.peek()));
dlg.show();
dlg.onClose.subscribe(() => barcode.dispose());
```

What the platform owns and what stays with u2:

- **Caption** — the platform's. `asDartInput` seeds it from the u2 label, then whatever binds the
  input (property name, `addCaption`) wins. `getInput()` hands over the u2 *editor*, not the
  input's root, because the Dart proxy wraps it in its own `ui-input-root` with a caption label.
- **Validation display** — the platform's. The u2 `validity` signal is registered as a platform
  validator, so the message shows once, in the platform's style; u2's own error line lives in the
  unmounted root and never appears.
- **`enabled` / `readOnly`** — set them on the object `asDartInput` returned, not on the u2 input
  (whose setter looks for controls under its own root, which no longer holds the editor).
- **Value** — u2's. The `value` signal stays the source of truth; read it with
  `input.value.peek()` or bind it. User edits fire `onInput` + `onChanged`; `setValue()` fires
  `onChanged` only, and nothing when the value is unchanged (echo suppression both ways).

## Lifecycle

`host()` makes a component a `DG.Widget` so the Dart kill-walk (view close) disposes its effects.
Pass `closeIn` when docking ad hoc — the pane ✕ fires `dockManager.onClosed` and never
`Widget.kill`:

```ts
const w = host(component, grok.shell.dockManager);
grok.shell.dockManager.dock(w.root, DG.DOCK_TYPE.RIGHT);
```

**Components that are never mounted are yours to dispose.** The kill-walk only reaches roots that
made it into the DOM — a lazily-built container (a collapsed accordion pane, a tab never opened)
never inserts the root, so nothing kills it. Same for value editors: the platform does not dispose
`DG.InputBase`s, so dispose the u2 input yourself when its host closes. Disposing the u2 input also
severs the bridge; `DartInput.detach()` severs only the bridge (it never owns the input).

## Pickers

```ts
const col = columnInput('Column', table, {filter: (c) => c.isNumerical});
const tbl = tableInput('Table');
```

Both are ordinary `ChoiceInput`s whose items follow the platform (`onColumnsChanged`,
`onTableAdded`/`onTableRemoved`) through the input's own scope — dispose the input and the
subscriptions go with it. The value is the **name**; the caller resolves it
(`table.col(col.value.peek())`, `grok.shell.table(tbl.value.peek())`), which keeps the value
serializable and the input identical to any other choice input.

## Schema-driven forms

`propertyForm(props, source)` builds a whole form out of `Property` metadata and edits `source` in
place. It is the primary surface: properties come from wherever the caller already has them —
`DG.Property.fromOptions`, a widget's or viewer's `getProperties()`, or the server
(`DomainRegistryClient.rowProperties()`).

```ts
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {propertyForm} from '@datagrok-libraries/u2/src/dg/index.js';

const props = [
  DG.Property.fromOptions({name: 'friendlyName', type: 'string', friendlyName: 'Name', nullable: false}),
  DG.Property.fromOptions({name: 'description', type: 'string', friendlyName: 'Description'}),
  DG.Property.fromOptions({name: 'language', type: 'string', friendlyName: 'Language',
    choices: ['python', 'r', 'javascript']}),
];

const script = await grok.dapi.scripts.find(id);
const form = propertyForm(props, script, {onChanged: (name) => console.log(`${name} edited`)});
form.addButtons((row) => row.append(ui.bigButton('Save', async () => {
  if (!form.validate())
    return;
  await grok.dapi.scripts.save(script);
  form.refresh();
})));
grok.shell.newView('Script', [form.root]);
```

**The form never saves anything.** `prop.set` mutates the object; persisting it is the caller's
call (`grok.dapi.*.save`, a domain `insert`), which keeps the platform's dapi flows — permissions,
validation, audit — exactly where they are. `refresh()` re-reads every field off the object
(echo-suppressed: no write-back, no `onChanged`) after a save or an external edit.

The generated form is inspectable: `form.properties` lists the rendered properties in layout
order, and `form.input(name)` returns the input editing a property — generated,
override-supplied, or a platform editor wrapped by `fromDartInput` — keyed by property name
regardless of its caption. From there, everything on `Input<T>` applies: `value`, `validity`,
`addValidator`, `enabled`.

How an input is picked, from `prop.propertyType` (or `type`):

| Metadata | Input |
|---|---|
| `choices` non-empty (any type) | `ChoiceInput` |
| `string`, `datetime` | `TextInput` — the date picker is v1's one deferred control, and an edit writes the typed string back, so override a `datetime` field until it lands |
| `int`, `bigint` | `NumberInput` int, `min`/`max` applied |
| `double`, `float`, `num`, `qnum` | `NumberInput` float, `min`/`max` applied |
| `bool` | `BoolInput` |
| anything else, or no `set` | disabled `TextInput` showing `String(value)` |

`caption ?? friendlyName ?? name` is the label, `name` is the input's name (so `form.getValues()`
is keyed by property name), `description` is the tooltip, and `nullable: false` adds a required
validator and drops the empty option from a choice. Everything else — `semType`, `units`,
`format`, `editor` — is ignored by the generator; that is what `overrides` is for:

```ts
const molecule = new TextInput({label: 'Structure'});   // or a real semType editor
const form = propertyForm(props, compound, {
  condensed: true,
  include: ['name', 'smiles', 'mw'],
  overrides: {
    smiles: {input: molecule},        // replaces generation entirely, still wired to prop.set
    mw: {label: 'MW, Da', inline: true},  // merges into the generated input's options
  },
});
```

A supplied `input` is laid out and wired (its value signal seeds from the property and writes
back), but the form never owns it — dispose it with whatever created it. Everything the form
generates it owns, so `form.dispose()` is the whole cleanup. A `bind` override works two-way; a
`value` override is ignored, since the object is the source of truth.

### `objectForm` and what can enumerate its own properties

`objectForm(source)` is `propertyForm(source.getProperties(), source)` — for objects that both
enumerate their properties and are the get/set target. That means **JS-implemented** `DG.Widget`
subclasses (a `JsViewer`, a custom filter or widget): their `addProperty` getters close over the
widget itself. A Dart-backed object — `DG.Viewer` and any other `DartWidget`, `DG.ViewBase` —
enumerates fine but reads and writes through its dart handle, so give `propertyForm` that handle
explicitly (it is what `ObjectPropertyBag` does internally):

```ts
const form = propertyForm(scatter.getProperties(), scatter.dart, {include: ['xColumnName', 'yColumnName']});
```

**`DG.Entity` is not one of them, and there is no generic way to enumerate an arbitrary entity's
properties from JS.** Both entity "get properties" calls return *values*, not metadata:
`Entity.getProperties()` and `grok.dapi.entities.getProperties(entity)` both resolve to a
name→value map of the entity's *dynamic* (schema-attached) properties (the latter's
`Map<Property, any>` return type in `dapi.ts` is stale — the Dart handler returns
`Map<String, dynamic>`). Every Dart `Entity` is a `PropMixin` and does have `properties`, but
js-api surfaces that (`grok_PropMixin_GetProperties`) only on widgets, views and viewers, and
the dynamic-property `Schema` has no JS data source at all. So an entity form declares its own
metadata — which is the honest thing anyway, since a form shows a chosen subset:

```ts
const values = await entity.getProperties();            // name -> value
const props = [
  DG.Property.fromOptions({name: 'Status', type: 'string', choices: ['New', 'Done']}),
  DG.Property.fromOptions({name: 'Priority', type: 'int', min: 1, max: 5}),
];
const form = propertyForm(props, values);
form.addButtons((row) => row.append(ui.bigButton('Save', () => entity.setProperties(values))));
```

Domain-schema rows are the one case where the server hands you real metadata —
`DomainRegistryClient.rowProperties('demo.issues')` returns `Property[]` with type, semType,
choices, min/max, nullable and the rendered `friendlyName`. Those getters are bound to the
**dart** row, so pass `DG.toDart(row)` as the source:

```ts
const props = await new DG.DomainRegistryClient().rowProperties('demo.issues');
const form = propertyForm(props, DG.toDart(row));
// ... form.validate(), then:
await grok.dapi.domains.table('demo.issues').insert(row.values);
```

which is the u2 counterpart of the platform's own `ui.input.form(toDart(row), properties)`.

### Relationship to func-param dialog generation

`objectForm` is the **JS-side** schema-driven generator and does not replace the Dart one. Function
parameter dialogs, `InputForm` and the property panel keep generating themselves in Dart from the
same `Property` metadata — they are wired into function calls, editors declared with `//editor:`,
and validation the platform runs server-side. A u2 control reaches those surfaces one input at a
time, through `asDartInput` and the `valueEditor` role (see above); replacing the dialog generator
itself is out of scope for u2. The two generators read the same metadata, so a property that
renders well in one renders well in the other.

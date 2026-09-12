# u2/dg — the platform layer

The only part of u2 that imports the platform surface of `datagrok-api` — eslint enforces it
(D3 in the u2 plan, `core/docs/features/ui2/PLAN.md`; since V-4 the core files re-export the
platform-free `datagrok-api/u2core`, the one allowed exception). Everything here bridges the
core onto Datagrok lifecycle, events, and value editors.

| Export | What it does |
|---|---|
| `toSignal` / `toObservable` | rxjs ↔ signals, both directions |
| `leakReport()` | live u2 scopes vs registered widgets |
| `appView({name, content, own?, ribbon?, toolbox?, status?})` | a u2 component tree as a platform view, riding the shell's own chrome; chrome controls and whatever `own` lists (the sources) are disposed with the content |
| `designerView(spec, options?)` | the spec designer as a view: live canvas, structure tree in the toolbox, Design/Run toggle, selection path in the status bar |
| `registerSpecNodeHandler()` | the `ObjectHandler` that renders a selected spec node (`SpecNodeRef`) in the context panel |
| `makeDesignerDroppable({element, active, onDragActive, onDrop})` | the platform's drag channel pointed at an element: a file or a function dragged out of Browse arrives as a `DropReading` (`readDrop`), each item ready for `dropNode` to make a `u2-func-source` of it |
| `asDartInput(input, {dataType?})` | u2 `Input` → `DG.JsInputBase` (a `DG.InputBase` everywhere) |
| `dartInputFor(build, {dataType?})` | the same, built from the bound property at `getInput()` time |
| `inputForProperty(prop, options?)` | the u2 editor a `Property` maps to, metadata and all |
| `columnInput(label, table, {filter?, rich?})` / `tableInput(label)` / `tablesInput(label)` | platform-bound choice inputs over columns and open tables |
| `columnRenderer(table)` | `ObjectRenderer<string>` over column names — type glyph, name, semType |
| `columnsInput(label, table, {filter?})` | many columns, picked in a searchable dialog |
| `columnsMapInput(label, keys, table)` / `aggregatedColumnsInput(label, table)` | key → column rows; column + aggregation rows |
| `metaInput(options?)` | read-only `DynamicInput` showing the object handler's card (Dart `MetaInput`) |
| `fileInput(label, {mode?, connection?})` / `filesInput(label, {accept?})` | one or many files, checked against the file share |
| `rsaInput(label, {accept?})` | a private key file, held masked |
| `propertyForm(props, source, options?)` / `objectForm(source, options?)` | a whole form out of `Property` metadata |
| `propertyEditor(record, options?)` | the editor of property metadata itself — the `IProperty` options |

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

## Property-aware value editors: `dartInputFor`

`asDartInput` wraps an input that already exists, so it can carry nothing the property says — a
registered `valueEditor` built that way loses every min/max/step/format/units/showSlider/
showPlusMinus the platform would have applied. `dartInputFor` defers construction instead:

```ts
grok.functions.register({
  signature: 'object u2IntValueEditor()',
  tags: 'valueEditor',
  options: {propertyType: 'int'},
  run: () => dartInputFor((prop) => inputForProperty(prop)),
});
```

**The timing contract.** `InputBase._forProperty` resolves the editor func, wraps the result in a
`JsInputProxy.fromFunc` and *immediately* runs `..caption = p.friendlyName ..bindProperty(p)`
(`input_base.dart:676-679`). Only afterwards does the proxy bind the JS object
(`jsi.dart = this`) and ask for its element from a `Timer` (`input_base_js_proxy.dart:18-25,27-38`).
So by the first `getInput()` the property is readable through `this.property`, and that is where
`build` runs — once. Consequences worth knowing:

- `getInput()` outside that path (a hand-built editor, a dialog) hands `build` a `null` property;
  `inputForProperty(null)` degrades to a text input.
- Values the platform writes before the editor exists (the proxy's `_tempValue` path) are buffered
  and applied to the input as programmatic writes.
- Repeated `getInput()` — the `FromJS` proxy and the `fromFunc` proxy each ask once — returns the
  same element and never rebuilds.
- The input is built under the bridge's own scope, so its effects are released by `detach()` (or by
  disposing the input, which severs the bridge) — a `Timer` callback has no ambient scope of its own.
- `u2` is `null` until the build, and so is the type it would name: `dataType` reads `object` until
  the editor exists. Afterwards the component id answers (`number-input` → `double`, `map-input` →
  `object`, …); only for a component outside that table does the type get inferred from the value —
  which reads `object` while the value is still null. Pass `dataType` explicitly in both cases.

`inputForProperty` is the mapping `propertyForm` uses, so a registered editor and a generated form
show the same control for the same metadata: bounds and step, `units` as the postfix, `format`
through `DG.format`, a clicker on bounded ints, a slider on floats or explicit `showSlider`, and a
`ChoiceInput` whenever the property carries `choices` — which matters because the valueEditor seam
is consulted *before* Dart's own choices branch (`input_base.dart:675` vs `:715-717`).

Pass `{assumeWritable: true}` from a value editor. A form renders a property it cannot write as a
disabled text box, and a `FuncParam` never carries a setter — so without it every func-param dialog
would get one. A value editor owns the value through the proxy, not through the property, so the
setter is beside the point there.

**Events and validation, both indirect.** The bridge fires `fireInput`/`fireChanged` from a
microtask rather than from the value effect: a platform subscriber that writes a u2 signal back —
the other half of a cross-synced form — would otherwise land inside the running effect round, and
the signal runtime aborts the round as a cycle and kills the effect for good. One edit is still one
pair of events, in that order, and a programmatic write (`setValue`) still fires `fireChanged`
alone; a value the input prunes on its own (`Input.system`, e.g. `ChoiceInput.setItems` dropping an
item that vanished) counts as programmatic too, as it does in Dart. Validators live on the *dart
handle*, and the platform re-homes the JS input on a second `JsInputProxy` when a `valueEditor` func
resolves — so the u2 validity feed is re-registered on every `getInput()`, once per handle. Without
that the form binds a proxy with no validators at all: the u2 editor shows its invalid state and the
dialog's OK stays enabled.

## Lifecycle

Since V-4 every `DG.Widget` IS a u2 `Control` (`DG.Widget extends DG.U2.Control`): the Dart
kill-walk (view close) disposes its scope, and disposing the scope kills the widget. To put a
plain component under the kill-walk, give its root the channel `appView` uses:

```ts
component.root.setAttribute('data-kill-on-close', 'true');
DG.Widget.registerCleanup(component.root, () => component.dispose());
grok.shell.dockManager.dock(component.root, DG.DOCK_TYPE.RIGHT);
```

Introspection is inherited, not delegated: `getProperties()`, `onEvent()`, `aiDescription` and
`getWidgetStatus()` are one surface on `DG.Widget`/`Control`, so the shell, the context panel and
the copilot interrogate a u2 component with the calls they already use on every platform widget.
(`host()` and its `grok.functions.register` minting were dropped in V-4.)

**Components that are never mounted are yours to dispose.** The kill-walk only reaches roots that
made it into the DOM — a lazily-built container (a collapsed accordion pane, a tab never opened)
never inserts the root, so nothing kills it. Same for value editors: the platform does not dispose
`DG.InputBase`s, so dispose the u2 input yourself when its host closes. Disposing the u2 input also
severs the bridge; `DartInput.detach()` severs only the bridge (it never owns the input).

## Pickers

```ts
const col = columnInput('Column', table, {filter: (c) => c.isNumerical});
const rich = columnInput('Column', table, {rich: true});   // searchable, with type glyphs
const tbl = tableInput('Table');
const tbls = tablesInput('Tables');                        // value: string[]
```

`columnInput`, `tableInput` and `tablesInput` are ordinary choice inputs whose items follow the
platform (`onColumnsChanged`, `onTableAdded`/`onTableRemoved`) through the input's own scope —
dispose the input and the subscriptions go with it. The value is the **name** (names, for the
multi-choice ones); the caller resolves it (`table.col(col.value.peek())`,
`grok.shell.table(tbl.value.peek())`), which keeps the value serializable and the input identical
to any other choice input. Whatever a picker holds when its column or table goes away is cleared.

`{rich: true}` swaps the native select for a type-ahead whose rows carry the column's type glyph,
name and semantic type — Dart's ColumnComboBox affordances, and a `ColumnPicker` rather than a
`ChoiceInput`, though the value contract is the same. The same rows come from `columnRenderer(table)`,
an `ObjectRenderer<string>` any u2 control can take.

`tableInput` and `tablesInput` carry the platform inputs' import action (`table_input.dart:37-58`,
`tables_multi_choice_input.dart:66-77`): a `folder-open` icon on the options rail opens a local
file, adds the table to the workspace and picks or checks it, as a user edit. It opens what js-api
can read on its own — csv/tsv/txt through `grok.data.parseCsv`, d42 through
`DG.DataFrame.fromByteArray` — and refuses anything else with the platform input's own
`File extension .x is not supported.`. Dart hands the bytes to `FileHandler`, which also opens
xlsx, sdf and whatever a package registered; js-api exposes no equivalent, so those formats stay
out of reach here.

## Column selectors

```ts
const cols = columnsInput('Features', table, {filter: (c) => c.isNumerical});  // value: string[]
const map = columnsMapInput('Map', ['id', {name: 'when', type: 'datetime'}], table);
const aggs = aggregatedColumnsInput('Aggregations', table);  // value: {column, aggregation}[]
```

`columnsInput` shows the summary Dart shows — `(3) age, sex, race`, or `(N) All` once everything
is in — and opens a searchable check-list with the picked columns on top; OK commits, cancel does
not, `changeTable(other)` starts over. The list is a u2 `Dialog` over a `VirtualList`, not the
platform's `selectColumns` modal, since the convergence retires that one; the additional-column
checks of the Dart input (`columns_input.dart:11-18,66-70`) are not ported and land when a consumer
(a pivot or aggregation UI migration) needs them.

`columnsMapInput` is a row per key — a plain column picker each, with the key's `type` restricting
what it offers (int and double count as one type, as in Dart). Unmapped keys stay out of the value.

`aggregatedColumnsInput` edits column + aggregation rows. The aggregation list follows the column's
type through `aggregationsFor(type)`, which is ddt's own per-type registry minus the exclusions
Dart's aggregation editor applies (`aggregated_columns_input.dart:8-12`) — so a string or bool
column offers the four counts, and a datetime one those plus `range`. Changing a row's column to a
type that does not offer the row's aggregation resets it to `defaultAggregation(column)`.

## Files and keys

```ts
const path = fileInput('Data', {mode: 'file'});                       // value: DG.FileInfo | null
const remote = fileInput('Data', {connection: myConnection});         // paths under one share
const batch = filesInput('Batch', {accept: ['.csv', '.sdf']});        // value: DG.FileInfo[]
const key = rsaInput('Private key');                                  // value: PEM text or BASE64
```

`fileInput` takes a path, a dropped file, or one chosen from disk. A typed path is checked against
the share — debounced 250 ms, one check at a time, a late answer to an abandoned path dropped —
and the input reports the platform's three states: empty (when not nullable), still checking, and
does not exist. A `connection` namespaces the path with its `nqName` before the check and drops
the local-file affordance, as the platform input does. What resolves comes out of the directory
listing, which is js-api's only per-path lookup: a path the server says exists but the listing does
not carry (the root of a connection, for one) stays valid and leaves the value as it was — js-api
cannot build a directory `FileInfo` the way `file_input.dart:147` does dart-side.

`filesInput` reads several files at once, each row carrying its own percentage and a ✕ that
cancels the read or drops the file; re-adding a name replaces it. The value — the loaded
`FileInfo`s — is published once every read has settled, failures included, and the input is
invalid while anything is loading or has failed.

`rsaInput` shows an upload button and a drop zone until a key is in and a masked field with a
Change button afterwards. The value is the file's text when it opens with `-----BEGIN`, and BASE64
of its bytes when it does not, so armoured and DER keys go through the same input. Keys dragged
out of the platform's file browser arrive as `FileInfo`s through `ui.makeDroppable` and must be
PEM; that registration is dart-side and has no undo, so it outlives `dispose()`.

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

How an input is picked — `inputType` first, then `editor`, then `propertyType` (or `type`), the
order `InputBase.forProperty` uses (`input_base.dart:668,689-714`):

| Metadata | Input |
|---|---|
| `inputType` | the editor it names, from the u2 half of the platform's factory map: `Text`, `Search`, `TextArea`, `Int`, `Float`, `BigInt`, `QNum`, `Slider`, `Bool`, `Switch`, `Date`, `Color`, `Font`, `Image`, `Radio`, `Choice`, `MultiChoice`, `Tags`, `List`, `Map`, `File`. An input type u2 has no editor for falls through to the rows below, where the platform would look for a JS-registered input |
| `editor` | `textarea` (or the `TextArea` spelling `Property.propertyOptions` itself uses) → `TextArea`, `password` → masked `TextInput`, `switch` → `BoolInput` switch, `slider` → `SliderInput` |
| `choices` non-empty (any scalar type) | `ChoiceInput` |
| `string` | `TextInput` |
| `datetime` | `DateTimeInput` |
| `int` | `NumberInput` int, `min`/`max` applied |
| `double`, `float`, `num` | `NumberInput` float, `min`/`max` applied |
| `bigint` | `BigIntInput` — the value is a JS `bigint`, which js-api marshals to and from a Dart `BigInt` (`wrappers_impl.ts:88,133`); a property that hands the digits over as text is read all the same |
| `qnum` | `QNumInput` — the value is the packed double the platform stores, `format` applied to its numeric part |
| `bool` | `BoolInput` |
| `list` | `MultiChoiceInput` with `choices`, `ListInput` (comma-separated) without |
| `map` | `MapInput` |
| `file`, `blob` | the dg `FileInput` — registered on `PlatformInputs` by `src/dg/index.ts`, since the router itself must load without the platform; importing only `object-form.js` leaves the file field read-only |
| anything else, or no `set` | disabled `TextInput` showing `String(value)` |

`caption ?? friendlyName ?? name` is the label, `name` is the input's name (so `form.getValues()`
is keyed by property name), `description` is the tooltip, and `nullable: false` adds a required
validator and drops the empty option from a choice. `semType` is the generator's own blind spot —
that is what `Editors` and `overrides` are for:

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
widget itself. A `DG.Viewer` enumerates fine but its properties are defined over its **look**
(`grok_Viewer_Get_Look(viewer.dart)`) — neither the viewer nor its dart handle is the get/set
target, and both throw — so use `viewerSettings(viewer)` (`dg/viewers/viewers.ts`), which hands
`propertyForm` the look and its user-editable properties:

```ts
const form = viewerSettings(scatter);   // propertyForm(userEditable props, grok_Viewer_Get_Look(scatter.dart))
```

The same fact drives the property tier: a widget's `propertyTarget` is what `prop.get`/
`prop.set` are handed — the look for a Dart viewer, `dart` for a `DartWidget`, the widget itself
otherwise (ordinary overrides in js-api `viewer.ts`/`widgets/base.ts` since V-4). The designer does not build its own editors for a viewer node:
`lookGrid(x, panel, funnel)` (`dg/designer/look-grid.ts`) mounts the platform `PropertyGrid` over
the live look with the viewer's frame (`grok_PropertyGrid_Update(pg, look, props, df.dart)` — the
frame is what gives `*ColumnName` rows a column picker) and turns each Design-mode commit
(`onPropertyValueChanged`) into one `set-prop` patch; an `object`-typed prop such as `filters` is
not captured — its write paths are *Add filter for column…*, a spec `set-prop`, or a live
`fg.props.filters = [...]`, which is the whole truth (panes are added, updated in place and
removed to match it). Note that `FilterGroup.props.columnNames` reads `[]` once panes exist —
the platform converts it into `filters`, so read `getOptions(true).look.filters`.

Other Dart-backed widgets (`DartWidget`, `DG.ViewBase`) read and write through their dart handle:
`propertyForm(widget.getProperties(), widget.dart)`.

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

## Editing the metadata itself: `propertyEditor`

`propertyEditor(record)` edits an `IProperty` record — what a property IS, not what it holds. Its
fields are the platform's own `Property.propertyOptions` vocabulary (js-api `property.ts:322-348`),
and their visibility follows the `applicableTo` the map declares: an identity section (name, type,
friendly name, description, nullable, semType, plus the `inputType`/`editor`/`units`/`format`/
`category` hints) that is always there, and a type-dependent one — `min`/`max`/`step`/`showSlider`/
`showPlusMinus` for a numerical property, `choices` for a string — rebuilt whenever the type or an
editor hint changes. The map itself is mirrored, not read: the control builds off-platform too.

The value model is a plain record, not a live property: `name` and `type` of a bound `DG.Property`
are not safely mutable, so **every edit reports the whole record** and the caller reconstructs.

```ts
const editor = propertyEditor(records['dose'], {onChanged: (options) => {
  records[options.name] = options;
  rebuild(DG.Property.fromOptions(options));   // and whatever hangs off it
}});
editor.setTarget(records['replicates']);       // re-target; the type-dependent fields rebuild
```

`setTarget` is the context-panel pattern: one editor, re-targeted per row, so nothing leaks per
click. `editor.options` is the same record as a signal, for a caller that would rather watch it.

Every field reports as it is typed except `name`, which reports on commit — blur or Enter. A caller
keys its records by that name (the example above does), and per keystroke `dose` → `name` would
rename the property through `n`, `na` and `nam` on the way. Pass `validators: {name: …}` to refuse a
name the caller cannot take: the message shows on the field, and the record keeps the name it had.

## Domain controls: `domains.table` and `u2.domain.*`

The platform adapter of the domain (EMS) seam (`src/dg/domain/`; recipe: `docs/recipes/crud-app.md`;
what both backends must agree on: `docs/domain-backend-contract.md`). `backends.domain` is filled
at dg import with `DgDomainBackend` over `grok.dapi.domains`, so a `DomainSource` (core) built
anywhere in the platform loads real rows; the same source over `MemoryDomainBackend` — a frame
host too, over `MemoryFrame` — is what the headless tests and the gallery run.

| Export | What it does |
|---|---|
| `domains.table(address)` | one await → `DomainTable`: `properties`, `info`, `access` (an `Access`), `handler` (the platform's `DomainObjectHandler` for the rows — never registered by u2), and the registries an app fills once: `actions` (`DomainAction {name, icon?, requires?, when?(row), run(row)}`), `validators.add(column, (value, row) => message \| null)`, `renderer` (the handler's rendering, the name column as caption) |
| `table.source(options)` / `table.draft(values)` | a started `DomainSource`; a `draft: true` source of one pristine draft. Both are tracked, so the controls find the handle (`DomainTable.of(src)`); `src.session` (core `SingleSession`) is the unit of work Save drives, `src.activate` the signal a list bumps on Enter |
| `domainForm(target, options?)` | `DomainForm` over a `DomainSource`, its `currentRow` signal, or a bare row (+ `source`): a fresh `propertyForm` per row under the row's access — the system columns left out (`system: 'footer'` shows them as a muted block: Id/Version/Created/Updated/Author, local times, the author resolved, a draft's id "assigned on save"), a readonly reference showing its target's name, a pristine draft's verdicts held back until a field is edited or left (`validate()` shows them all), an edited field marked with an amber edge, `focus()` on the first field; validity = schema rules + `table.validators` + the writer's per-cell errors (over the platform the editor maps the server's column errors onto the cells and resolves a 409 itself). Paired through the source: it guards `src.save()` (a refusal is "Cannot save: <Caption> is required" whichever path ran — the button, Ctrl+S, `cmd:save`), takes the focus back on the session's `onSaved`/`onDiscarded` and on `src.activate`; a refused save is one balloon, a 403 also drops the access caches |
| `domainList(src, options?)` | `DomainList`: a `VirtualList` over `src.rows`, `brief` or `cards` (the recipe card: title, one muted description line, the creation time; a handler's own card where it defines one), Open (saved rows) + Delete (`requires: 'delete'`) + the table's and the list's own actions per row; a row marked deleted stays struck through with Restore until the save; selection = `src.currentRow` (a loaded list starts on its first row), Enter bumps `src.activate` (the paired form takes the focus), Delete deletes, the selected row's actions rove into the tab order; the next page near the bottom, loading/empty/error under the rows |
| `domainPick(table, options?)` | `DomainPick`, an `Input<string \| null>` face over a `TypeAhead` that queries the target table by its name column (else the first business-key column, else the id); `PickInput` is the same face over any type-ahead — the `User`/`Group` editors ride it |
| `saveButton(src \| session)` / `discardButton(src \| session)` / `newButton(src, values?)` | buttons over a `DomainSession` (a source stands for its own): disabled while clean or saving; the session balloons "<Singular> saved" once and emits `onSaved`, so no form is passed around; New adds a pristine draft (`src.newRow`), hidden without `insert` |
| `registerDomainComponents(reg)` | the `u2-domain-form` / `-list` / `-pick` tags (also run by `registerPlatformComponents`); a control binds to its source whole: `bind: {source: "$.issues.source"}` |

Editor rules registered at import: a property whose `semType` is a table address
(`<schema>.<table>`) gets `domainPick` over that table; `User` → the user picker, `Group` → the
group picker — everywhere `propertyForm` is generated, not only in domain forms.

### Access

`Access` (core) is the one protocol every control consults; the server answers it in exactly that
shape (`DomainTableClient.access()` → `{can, fields}`) and nothing asks the server twice.

- `can(capability)`: `view`, `insert`, `edit`, `delete`, `share` (plus a schema's custom names).
  `Action.requires` names one; a denied action is not rendered.
- `field(name)`: `hidden` where the server did not list the column, `editable` where it listed
  it so AND the row may be written, `readonly` otherwise. The write capability is `edit` for an
  existing row and `insert` for a draft (`access.forDraft()`) — a listed `editable` only means
  "not column-restricted", never "this caller may write this row".
- `row(r)`: the per-row view. A `withAccess` query adds `~can_edit`, `~can_delete`, `~can_share`
  to every row (`Access.ROW_COLUMNS`, the documented server set) — computed by the same predicates
  the server enforces on write; `row()` reads any `~can_<name>` the table has a capability for. A
  boolean the row carries REPLACES the table's answer for that capability either way: a row-mode
  table's table-level flags are false negatives, so a row may be editable under a table that says
  `edit: false`, and read-only under one that says `true`. A column that is absent or not a
  boolean is not carried — off row mode the server sends `~can_share` as null and omits it from
  a frame. `insert` is never per row.
- The controls pick the view for you: `domainForm` renders a draft under `forDraft()` and an
  existing row under `row(r)`; `domainList` filters each row's actions under `row(r)` (a draft
  under the table's access — its `~can_*` cells are the frame's defaults, not truth). The frame
  path hides the three columns from every export the way the editor hides its own service columns.

### What the frame host honours (STATE-CONTRACT)

Over the platform backend a source is a frame with the js-api `DomainFrameEditor` attached as its
single writer (`EditorEditState`); over the memory backend a `MemoryFrame` with `MemoryEditState`.
Rows are keyed the way `FrameRows` keys them — the `id` cell, `Rows.draftKey(index)` (`~row:<index>`)
for a draft that has no id yet; `Rows` (core) is the one holder of the `~` conventions (`Rows.STATE`,
`isService`, `isDraft`). What that means for an app:

- **H6** — a programmatic re-query (a bound `query` changing) is skipped while the source is dirty;
  a user-initiated one must go through the unsaved-changes gate. A refresh drops the batch.
- **H7** — never export or hand off the source's frame: it carries the editing state and the
  access columns (export-tagged, but that is the safety net). Query a fresh frame.
- **H8** — a re-query replaces the frame and detaches the editor from the old one; the next page
  (`loadMore`) is appended into the same frame, pending edits kept.
- **H10** — everything mirrored (`isDirty`, `changeCount`, `validity`, `errorOf`) follows the
  editor's `onChanged`, never frame events; a direct `df.set` is invisible to save.

/* The Design step's right pane: one panel per selection kind, rebuilt when the selection or the
   model changes as keyed parts (title, form, sections), each in its own scope: a part that comes out
   of the rebuild as it is shown — the same make-up, its controls holding what the user left in them —
   stays, so the focus and the caret in it survive; a changed part is replaced, and focus inside it is
   put back on the element with the same keys. Schema:
   name, friendly name, writable, the access grid for every table. Table: logical and friendly
   names, the key (read-only), the name and searchable pickers, the read-only opt-out under a
   writable binding, the relationships with an "include <target>" link where including the
   target would make a ref, the access grid with the schema's rows inherited. Column: logical
   name, the mapped type beside the warehouse's, required, name, searchable, its reference with
   its status, visibility. The field offer follows the editor context: `view` renders every
   field as text, never as a dead input. Access is modelled beside the manifest and rendered
   here; the dialog applies it after Create. */
import {Control, Scope, signal, ReadonlySignal, div, link, span, each, focusMark, refocus, Form, Section, TextInput,
  BoolInput, ChoiceInput, RadioInput, ChipsInput, TagsInput, NumberInput, AccessGrid, badge, camelCaseToWords}
  from '@datagrok-libraries/u2';
import type {Input, AccessRow, InheritedAccessRow} from '@datagrok-libraries/u2';
import {ObjectForm} from '@datagrok-libraries/u2/src/dg/forms/object-form.js';
import type {EditorContext, FieldOffer} from './editor-context.js';
import {fieldOffer} from './editor-context.js';
import {AccessModel, ManifestModel, EXTEND_TYPES} from './manifest-model.js';
import type {AccessGrant, AccessPrincipal, AccessScope, ColumnAttribute, ColumnView, DriftView, ManifestDiagnostic,
  ManifestSelection, RelationView, TableView} from './manifest-model.js';
import {ManifestRules} from './manifest-rules.js';
import {ManifestTree} from './manifest-tree.js';

/** Builds a look-up control that hands every picked principal to `onPick` — the platform's group
 * and user picker, in the access grids and beside the visibility chips. */
export type PrincipalPicker = (onPick: (principal: AccessPrincipal) => void) => HTMLElement;

export interface ManifestContextPanelOptions {
  selected: ReadonlySignal<ManifestSelection>;
  access: AccessModel;
  context: EditorContext;
  /** The groups the access pickers offer up front. */
  groups?: AccessPrincipal[];
  /** The look-up that finds any other; without it the pickers offer `groups` alone. */
  principalPicker?: PrincipalPicker;
  /** Diagnostics addressed by manifest path; the ones about the selected node are listed on top. */
  diagnostics?: ReadonlySignal<ManifestDiagnostic[]>;
  /** Why the Writable switch cannot be turned on here; the switch is then disabled with this hint. */
  writableDisabled?: string;
  /** Selects a node in the tree — where a column the panel added is taken to. */
  onSelect?: (selection: ManifestSelection) => void;
}

const SCHEMA: AccessScope = {kind: 'schema'};

const CAPABILITIES = [{name: 'view', label: 'View'}, {name: 'edit', label: 'Edit'}, {name: 'delete', label: 'Delete'}];
const CREATOR = 'You (creator)';
const EVERYONE = 'Everyone who may see the row';
const SOME = 'Only these groups';
const NONE = '— none —';

/** Attributes that point at minted ids or group names: not part of what a part shows. */
const REFS = ['id', 'for', 'name', 'aria-labelledby', 'aria-describedby', 'aria-controls', 'aria-owns',
  'data-u2-owner-id'];

/** A keyed piece of the panel, built in its own scope. */
interface Part {
  key: string;
  el: HTMLElement;
  scope: Scope;
  /** Its make-up as built, values aside. */
  shape: string;
  /** What its controls held as built. */
  values: string;
}

export class ManifestContextPanel extends Control {
  readonly offer: FieldOffer;

  private readonly _access: AccessModel;
  private readonly _groups: AccessPrincipal[];
  private readonly _picker: PrincipalPicker | undefined;
  private readonly _diagnostics: ReadonlySignal<ManifestDiagnostic[]> | undefined;
  private readonly _writableDisabled: string | undefined;
  private readonly _onAdded: ((selection: ManifestSelection) => void) | undefined;
  private readonly _shapes = new WeakMap<HTMLElement, string>();

  constructor(readonly model: ManifestModel, options: ManifestContextPanelOptions) {
    super();
    this.offer = fieldOffer(options.context);
    this._access = options.access;
    this._groups = options.groups ?? [];
    this._picker = options.principalPicker;
    this._diagnostics = options.diagnostics;
    this._writableDisabled = options.writableDisabled;
    this._onAdded = options.onSelect;
    this.root.classList.add('u2-manifest-panel');
    this.root.dataset.u2 = 'manifest-panel';
    this.runInScope(() => each(this.root, () => this._parts(options.selected.value), (p) => p.key, (p) => {
      Scope.ambient!.own(() => p.scope.dispose());
      this._shapes.set(p.el, p.shape);
      return p.el;
    }, {update: (el, p) => {
      if (this._shapes.get(el) !== p.shape || ManifestContextPanel._values(el) !== p.values)
        return false;
      p.scope.dispose();
      return true;
    }}));
  }

  /** The parts the selection shows now, keyed by the selection and the part. Focus on an element a
   * replaced part takes away goes to the element with the same keys in its replacement. */
  private _parts(selection: ManifestSelection): Part[] {
    this.model.revision.value;
    this.model.name.value;
    const active = document.activeElement;
    const mark = focusMark(this.root);
    if (mark !== null) {
      queueMicrotask(() => {
        if (!this.scope.isDisposed && active !== null && !active.isConnected)
          refocus(this.root, mark);
      });
    }
    const problems = (this._diagnostics?.value ?? [])
      .filter((d) => ManifestModel.sameSelection(this.model.resolvePath(d.path), selection));
    const parts = selection.kind === 'schema' ? this._schema() :
      selection.kind === 'table' ? this._table(selection.table) :
        this._column(selection.table, selection.column);
    if (problems.length > 0) {
      parts.splice(1, 0, this._part('diagnostics', () => div(problems.map((p) =>
        div([badge(p.code, {variant: 'error'}), span(p.message)], 'u2-manifest-panel-diagnostic')),
      'u2-manifest-panel-diagnostics')));
    }
    const at = selection.kind === 'schema' ? 'schema' : selection.kind === 'table' ? `table:${selection.table}` :
      `column:${selection.table}.${selection.column}`;
    for (const p of parts)
      p.key = `${at}/${p.key}`;
    return parts;
  }

  private _part(key: string, build: () => HTMLElement | Control): Part {
    const scope = new Scope(this.scope, {label: 'part', lazy: true});
    const built = Scope.runWith(scope, build);
    const el = Control.is(built) ? built.root : built;
    return {key, el, scope, shape: ManifestContextPanel._shape(el), values: ManifestContextPanel._values(el)};
  }

  private static _shape(el: Element): string {
    const attrs = el.getAttributeNames().filter((n) => !REFS.includes(n)).map((n) => `${n}="${el.getAttribute(n)}"`);
    const type = el.tagName === 'INPUT' ? ` type="${(el as HTMLInputElement).type}"` : '';
    let shape = `<${el.tagName}${type} ${attrs.join(' ')}>`;
    for (const child of el.childNodes) {
      shape += child.nodeType === 1 ? ManifestContextPanel._shape(child as Element) :
        child.nodeType === 3 ? child.textContent : '';
    }
    return shape;
  }

  private static _values(el: Element): string {
    return JSON.stringify([...el.querySelectorAll<HTMLInputElement>('input, select, textarea')]
      .map((c) => [c.value, c.checked === true]));
  }

  private _schema(): Part[] {
    const model = this.model;
    const name = model.name.peek();
    const writable = model.writable.peek();
    if (model.extending) {
      return [this._part('title', () => ManifestContextPanel._title('Schema', name)),
        this._part('note', () => ManifestContextPanel._note('Pick a table to add properties of your own to it. ' +
          'Columns with a lock come from the plugin and are read-only.'))];
    }
    return [this._part('title', () => ManifestContextPanel._title('Schema', name)),
      this._part('form', () => this._schemaForm(name, writable)),
      this._part('access', () => this._schemaAccess(writable))];
  }

  private _schemaForm(name: string, writable: boolean): Form {
    const model = this.model;
    const form = new Form({layout: 'wide'});
    const friendly = model.friendlyName.peek();
    this._field(form, 'Friendly name', 'friendlyName', friendly, () => new TextInput({label: 'Friendly name',
      name: 'friendlyName', value: friendly, commitOn: 'change', onChanged: (v) => model.setSchemaFriendlyName(v)}));
    this._field(form, 'Identifier', 'name', name, (hint) => {
      const input = new TextInput({label: 'Identifier', name: 'name', value: name, commitOn: 'change',
        onChanged: (v) => model.setSchemaName(v), ...hint});
      input.addValidator((v) => model.checkSchemaName(v));
      return input;
    }, `registered as ext_${name}`, this.offer.identifier);
    if (this.offer.description) {
      const description = model.description.peek();
      this._field(form, 'Description', 'description', description, () => new TextInput({label: 'Description',
        name: 'description', value: description, commitOn: 'change', onChanged: (v) => model.setDescription(v)}));
    }
    if (this.offer.writable) {
      const hint = this._writableDisabled;
      this._field(form, 'Writable', 'writable', writable ? 'Yes' : 'No', (extra) => new BoolInput({label: 'Writable',
        name: 'writable', value: writable, enabled: hint === undefined, tooltipText: hint,
        onChanged: (v) => model.setWritable(v), ...extra}),
      hint ?? 'users with Edit on a table may insert, update and delete rows in the warehouse');
    }
    form.addElement(ManifestContextPanel._note('Queries run as the platform service with the connection\'s ' +
      'stored credentials; users need View on a table, nothing on the connection.'));
    if (model.editing && model.catalog === 'unknown') {
      form.addElement(ManifestContextPanel._note('The warehouse catalog could not be read: the registered tables ' +
        'and columns are shown as they were, nothing is marked missing, and no new table is offered.'));
    }
    return form;
  }

  private _schemaAccess(writable: boolean): Section {
    const model = this.model;
    const access = new Section({title: 'Access — every table', collapsible: false});
    const locked = writable ? [] : ['edit', 'delete'];
    if (this._access.editing) {
      const tables = model.tables.peek().filter((t) => t.included).map((t) => t.remote);
      access.add(this._accessGrid('schema', this._access.everyTable(tables), [], locked,
        (rows) => this._access.setEveryTable(tables, rows)));
      access.add(ManifestContextPanel._note('What every included table you may share grants alike; a change here ' +
        'is written on each of them, and a table\'s own rows are edited on its panel.'));
    } else {
      access.add(this._accessGrid('schema', this._access.grantsOf(SCHEMA).peek(), [ManifestContextPanel._creator()],
        locked, (rows) => this._access.setGrants(SCHEMA, rows)));
      access.add(ManifestContextPanel._note('Applied after Create as the same grant on every included table — ' +
        'there is no schema-wide row access. Edit and Delete need a writable binding.'));
    }
    return access;
  }

  private _table(remote: string): Part[] {
    const table = this.model.table(remote);
    const title = this._part('title', () => {
      const title = ManifestContextPanel._title('Table', remote);
      if (table === undefined)
        return title;
      if (!table.included) {
        const status = !table.bindable ? ManifestTree.shortReason(table.code) :
          table.drafted ? (table.registered ? 'removed' : 'not included') : 'not in this draft';
        title.append(badge(status, {variant: table.bindable ? 'warning' : 'error'}),
          span(!table.bindable || !table.drafted ? table.reason ?? '' :
            table.registered ? 'unregistered on Save — its grants, restrictions and saved filters go with it' :
              'check it in the tree to expose it', 'u2-manifest-node-hint'));
      } else
        ManifestContextPanel._drift(title, table.drift);
      return title;
    });
    if (table === undefined || !table.included && (!table.bindable || !table.drafted))
      return [title];
    if (this.model.extending)
      return [title, this._part('form', () => this._extendTableForm(remote, table))];
    // a registered table taken out is on its way to the purge: its fields are read, not edited
    const live = table.included || !table.registered;
    return [title, this._part('form', () => this._tableForm(remote, table, live)),
      this._part('relations', () => this._tableRelations(remote)),
      this._part('access', () => this._tableAccess(remote, table, live))];
  }

  private _tableForm(remote: string, table: TableView, live: boolean): Form {
    const model = this.model;
    const columns = model.columns(remote).peek();
    const form = new Form({layout: 'wide'});
    this._field(form, 'Logical name', 'logical', table.logical, () => {
      const input = new TextInput({label: 'Logical name', name: 'logical', value: table.logical, commitOn: 'change',
        onChanged: (v) => model.renameTable(remote, v)});
      input.addValidator((v) => model.checkTableName(remote, v));
      return input;
    }, table.registered ? 'registered — a logical name is for life' : undefined, !table.registered);
    // an unsaid friendly name is the derived one: shown, not sent
    const derived = camelCaseToWords(table.logical);
    this._field(form, 'Friendly name', 'friendlyName', table.friendlyName || derived, () => new TextInput({
      label: 'Friendly name', name: 'friendlyName', value: table.friendlyName, placeholder: derived, commitOn: 'change',
      onChanged: (v) => model.setFriendlyName(remote, v)}), undefined, live);
    const key = table.key.map((k) => columns.find((c) => c.remote === k)?.logical ?? k).join(', ');
    const keyField = ObjectForm.readonlyField('Key', 'key', key);
    keyField.value.append(ManifestContextPanel._hint('the primary key; the row id encodes it'));
    form.addElement(keyField.row);
    const strings = columns.filter((c) => c.supported && c.included && c.type === 'string')
      .map((c) => ({value: c.remote, label: c.logical}));
    const pick = (label: string, name: string, current: ColumnView | undefined,
      apply: (remote: string | null) => void): void => {
      this._field(form, label, name, current?.logical ?? NONE, () => new ChoiceInput({label, name, items: strings,
        value: current?.remote ?? null, onChanged: apply}), undefined, live);
    };
    if (this.offer.nameColumn)
      pick('Name column', 'nameColumn', columns.find((c) => c.isName), (r) => model.setNameColumn(remote, r));
    if (this.offer.searchable)
      pick('Searchable', 'searchable', columns.find((c) => c.searchable), (r) => model.setSearchable(remote, r));
    if (this.offer.writable && model.writable.peek()) {
      this._field(form, 'Read-only', 'readOnly', table.readOnly ? 'Yes' : 'No', (hint) => new BoolInput({
        label: 'Read-only', name: 'readOnly', value: table.readOnly, onChanged: (v) => model.setReadOnly(remote, v),
        ...hint}),
      'opt this table out of writes', live);
    }
    return form;
  }

  private _tableRelations(remote: string): Section {
    const model = this.model;
    const relations = new Section({title: 'Relationships', collapsible: false});
    const rels = model.relations.peek().filter((r) => r.table === remote);
    if (rels.length === 0) {
      const unread = (this._diagnostics?.value ?? []).some((d) => d.code === 'external-relations-unavailable');
      relations.add(ManifestContextPanel._note(unread ?
        'The foreign keys of this schema could not be read — every column stays a plain value.' :
        'The warehouse reports no foreign keys on this table.'));
    }
    for (const r of rels)
      relations.add(this._relation(r, `${r.column} → ${r.targetTable}`));
    return relations;
  }

  private _tableAccess(remote: string, table: TableView, live: boolean): Section {
    const model = this.model;
    const scope: AccessScope = {kind: 'table', table: remote};
    const rows = this._access.grantsOf(scope).peek();
    const locked = model.writes(table) ? [] : ['edit', 'delete'];
    const access = new Section({title: 'Access — this table', collapsible: false});
    if (this._access.editing) {
      const editable = this._access.canEdit(scope);
      access.add(this._accessGrid(remote, rows, [], locked, (changed) => this._access.setGrants(scope, changed),
        {enabled: editable && live}));
      if (!editable) {
        access.add(ManifestContextPanel._note(table.registered ?
          'You cannot share this table: its grants cannot be read here and stay as they are.' :
          'Access cannot be read — it stays as it is.'));
      }
      for (const g of rows.filter((g) => (g.other ?? []).length > 0))
        access.add(ManifestContextPanel._note(`${g.group.label} also holds ${g.other!.join(', ')} — kept as is.`));
    } else {
      const schemaGrants = this._access.grantsOf(SCHEMA).peek();
      const schemaRows = schemaGrants
        .map((g): InheritedAccessRow => ({...ManifestContextPanel._row(g), from: '(every table)'}));
      access.add(this._accessGrid(remote, rows, [ManifestContextPanel._creator(), ...schemaRows], locked,
        (changed) => this._access.setGrants(scope, changed), {inheritedGroups: schemaGrants.map((g) => g.group)}));
    }
    return access;
  }

  private _column(table: string, remote: string): Part[] {
    const column = this.model.column(table, remote);
    const title = this._part('title', () => {
      const title = ManifestContextPanel._title('Column', `${table}.${remote}`);
      if (column !== undefined && !column.supported)
        title.append(badge('not bindable', {variant: 'error'}), span(column.reason ?? '', 'u2-manifest-node-hint'));
      else if (column !== undefined)
        ManifestContextPanel._drift(title, column.drift);
      return title;
    });
    if (column === undefined || !column.supported)
      return [title];
    if (this.model.extending) {
      return [this._part('title', () => this._extendColumnTitle(table, remote, column)),
        this._part('form', () => column.owned && column.included ? this._extendColumnForm(table, remote, column) :
          this._pluginColumnForm(column))];
    }
    const parts = [title, this._part('form', () => this._columnForm(table, remote, column))];
    if (column.relations.length > 0) {
      parts.push(this._part('reference', () => {
        const reference = new Section({title: 'Reference', collapsible: false});
        for (const r of column.relations)
          reference.add(this._relation(r, `foreign key to ${r.targetTable}`));
        return reference;
      }));
    }
    if (this.offer.access)
      parts.push(this._part('visibility', () => this._visibility(table, remote, column)));
    return parts;
  }

  // ─────────────────────── extend mode ───────────────────────

  private _extendTableForm(remote: string, table: TableView): Form {
    const columns = this.model.columns(remote).peek();
    const own = columns.filter((c) => c.owned && c.included).length;
    const form = new Form({layout: 'wide'});
    form.addElement(ObjectForm.readonlyField('Friendly name', 'friendlyName',
      table.friendlyName || camelCaseToWords(table.logical)).row);
    form.addElement(ObjectForm.readonlyField('Properties', 'properties',
      `${columns.length - own} from the plugin, ${own} of your own`).row);
    if (this.offer.editable) {
      const add = link('+ Add property', () => this._add(remote));
      add.dataset.u2 = 'manifest-add-column';
      form.addElement(add);
    }
    return form;
  }

  /** The column's title: who owns it, what the save does to it, and the way to remove or restore it. */
  private _extendColumnTitle(table: string, remote: string, column: ColumnView): HTMLElement {
    const title = ManifestContextPanel._title('Column', `${table}.${column.logical}`);
    if (!column.owned) {
      title.append(badge('plugin'));
      return title;
    }
    if (!column.registered)
      title.append(badge('new', {variant: 'accent'}));
    else if (!column.included)
      title.append(badge('removed', {variant: 'error'}), span('dropped on Save, with its values', 'u2-manifest-node-hint'));
    else
      title.append(badge('custom', {variant: 'success'}));
    if (!this.offer.editable)
      return title;
    let action: HTMLElement;
    if (column.included) {
      action = link('Remove', () => this.model.removeColumn(table, remote));
      action.dataset.u2 = 'manifest-remove-column';
    } else {
      action = link('Restore', () => this.model.includeColumn(table, remote, true));
      action.dataset.u2 = 'manifest-restore-column';
    }
    action.classList.add('u2-manifest-panel-remove');
    title.append(action);
    return title;
  }

  /** A plugin column, or one of the user's own on its way out: every field as text. */
  private _pluginColumnForm(column: ColumnView): Form {
    const a = column.attributes;
    const form = new Form({layout: 'wide'});
    const text = (label: string, name: string, value: string | undefined): void => {
      if (value !== undefined && value !== '')
        form.addElement(ObjectForm.readonlyField(label, name, value).row);
    };
    text('Name', 'logical', column.logical);
    text('Friendly name', 'friendlyName', a.friendlyName);
    text('Type', 'type', column.type === 'ref' ? `ref → ${column.relations.find((r) => r.ref)?.targetLogical ?? ''}` :
      column.type);
    text('Description', 'description', a.description);
    text('Choices', 'choices', a.choices?.join(', '));
    text('Range', 'range', a.min === undefined && a.max === undefined ? undefined : `${a.min ?? ''} – ${a.max ?? ''}`);
    text('Format', 'format', a.format);
    text('Default', 'default', a.default === undefined ? undefined : String(a.default));
    if (!column.owned) {
      form.addElement(ManifestContextPanel._note('Defined by the plugin, so it cannot be changed here; plugin ' +
        'upgrades keep it in step.'));
    }
    return form;
  }

  /** A column of the user's own: the declaration's descriptive keys, by type; the name and the type
   * only until it is saved. */
  private _extendColumnForm(table: string, remote: string, column: ColumnView): Form {
    const model = this.model;
    const a = column.attributes;
    const set = (key: ColumnAttribute) => (v: unknown): void => model.setColumnAttribute(table, remote, key, v);
    const form = new Form({layout: 'wide'});
    form.add(new TextInput({label: 'Friendly name', name: 'friendlyName', value: a.friendlyName ?? '',
      placeholder: camelCaseToWords(column.logical), commitOn: 'change', onChanged: set('friendlyName')}));
    this._field(form, 'Name', 'logical', column.logical, (hint) => {
      const input = new TextInput({label: 'Name', name: 'logical', value: column.logical, commitOn: 'change',
        onChanged: (v) => model.renameColumn(table, remote, v), ...hint});
      input.addValidator((v) => model.checkColumnName(table, remote, v));
      return input;
    }, column.registered ? 'saved — used in filters and links' : 'fixed once saved', !column.registered);
    this._field(form, 'Type', 'type', column.type, (hint) => new ChoiceInput({label: 'Type', name: 'type',
      items: [...EXTEND_TYPES], value: column.type, nullable: false,
      onChanged: (v) => model.setColumnType(table, remote, v as string), ...hint}),
    column.registered ? 'saved — a new type would rewrite the stored values' : 'fixed once saved', !column.registered);
    form.add(new TextInput({label: 'Description', name: 'description', value: a.description ?? '', commitOn: 'change',
      onChanged: set('description')}));
    const type = column.type;
    if (ManifestRules.attributeFits('choices', type)) {
      form.add(new TagsInput<string>({label: 'Choices', name: 'choices', value: a.choices ?? [], allowNew: true,
        placeholder: 'Type a value, Enter to add', onChanged: set('choices')}));
    }
    if (ManifestRules.attributeFits('min', type)) {
      const mode = type === 'int' ? 'int' : 'float';
      form.add(new NumberInput({label: 'Min', name: 'min', value: a.min ?? null, mode, onChanged: set('min')}));
      form.add(new NumberInput({label: 'Max', name: 'max', value: a.max ?? null, mode, onChanged: set('max')}));
    }
    if (ManifestRules.attributeFits('format', type)) {
      form.add(new TextInput({label: 'Format', name: 'format', value: a.format ?? '', commitOn: 'change',
        placeholder: type === 'datetime' ? 'yyyy-MM-dd' : '0.00', onChanged: set('format')}));
    }
    if (type === 'string' && (a.choices ?? []).length > 0) {
      form.add(new ChoiceInput({label: 'Default', name: 'default', items: a.choices!,
        value: typeof a.default === 'string' && a.choices!.includes(a.default) ? a.default : null,
        onChanged: set('default')}));
    } else if (type === 'string')
      form.add(new TextInput({label: 'Default', name: 'default', value: String(a.default ?? ''), commitOn: 'change', onChanged: set('default')}));
    else if (type === 'int' || type === 'float') {
      form.add(new NumberInput({label: 'Default', name: 'default', value: typeof a.default === 'number' ? a.default : null,
        mode: type, onChanged: set('default')}));
    } else if (type === 'bool') {
      form.add(new ChoiceInput({label: 'Default', name: 'default', items: ['true', 'false'],
        value: a.default === undefined ? null : String(a.default),
        onChanged: (v) => set('default')(v === null ? null : v === 'true')}));
    }
    form.add(new BoolInput({label: 'Required', name: 'required', value: false, enabled: false,
      postfix: 'properties of your own stay optional, so rows saved before stay valid'}));
    return form;
  }

  private _add(table: string): void {
    const column = this.model.addColumn(table);
    if (column !== null)
      this._onAdded?.({kind: 'column', table, column});
  }

  private _columnForm(table: string, remote: string, column: ColumnView): Form {
    const model = this.model;
    const form = new Form({layout: 'wide'});
    this._field(form, 'Logical name', 'logical', column.logical, () => {
      const input = new TextInput({label: 'Logical name', name: 'logical', value: column.logical, commitOn: 'change',
        onChanged: (v) => model.renameColumn(table, remote, v)});
      input.addValidator((v) => model.checkColumnName(table, remote, v));
      return input;
    }, column.registered ? 'registered — a logical name is for life' : undefined, !column.registered);
    const type = column.type === 'ref' ? `ref → ${column.relations.find((r) => r.ref)?.targetLogical ?? ''}` :
      column.type;
    const typeField = ObjectForm.readonlyField('Type', 'type',
      column.dbType === undefined ? type : `${type} ← ${column.dbType}`);
    if (column.drift !== undefined && column.drift.kind !== 'unknown')
      typeField.value.append(ManifestContextPanel._hint(column.drift.reason));
    form.addElement(typeField.row);
    if (this.offer.required) {
      this._field(form, 'Required', 'required', column.required ? 'Yes' : 'No', (hint) => new BoolInput({
        label: 'Required', name: 'required', value: column.required, enabled: !column.isKey,
        onChanged: (v) => model.setRequired(table, remote, v), ...hint}),
      column.isKey ? 'key columns are always required' :
        'NOT NULL is not reported by the warehouse; tick it if you know');
    }
    const isString = column.type === 'string';
    if (this.offer.nameColumn) {
      this._field(form, 'Name column', 'isName', column.isName ? 'Yes' : 'No', (hint) => new BoolInput({
        label: 'Name column', name: 'isName', value: column.isName, enabled: isString,
        onChanged: (v) => model.setNameColumn(table, v ? remote : null), ...hint}),
      'shown wherever a row is referred to');
    }
    if (this.offer.searchable) {
      this._field(form, 'Searchable', 'searchable', column.searchable ? 'Yes' : 'No', (hint) => new BoolInput({
        label: 'Searchable', name: 'searchable', value: column.searchable, enabled: isString,
        onChanged: (v) => model.setSearchable(table, v ? remote : null), ...hint}), 'one per table');
    }
    return form;
  }

  private _visibility(table: string, remote: string, column: ColumnView): Section {
    const section = new Section({title: 'Visibility', collapsible: false});
    if (column.isKey) {
      section.add(ManifestContextPanel._note(
        'A key column is visible to everyone who may see the row: the row id encodes it.'));
      return section;
    }
    const access = this._access;
    const entry = access.columnOf(table, remote);
    const current = entry?.groups ?? null;
    if (!this.offer.editable || !access.canEditColumn(table, remote) || entry?.unknown === true) {
      section.add(ManifestContextPanel._note(entry?.unknown === true ?
        `${SOME} — which ones cannot be read here; it stays as it is` :
        current === null ? EVERYONE :
          `${SOME}: ${current.length === 0 ? 'none yet' : current.map((g) => g.label).join(', ')}`));
      if (this.offer.editable && !access.canEditColumn(table, remote))
        section.add(ManifestContextPanel._note('You cannot share this table\'s columns: the visibility stays as it is.'));
      return section;
    }
    if (entry?.edit !== undefined && entry.edit.length > 0)
      section.add(ManifestContextPanel._note(`May also edit it: ${entry.edit.map((g) => g.label).join(', ')}.`));
    const some = signal(current !== null);
    const radio = new RadioInput({label: 'Visible to', name: 'visibility', items: [EVERYONE, SOME],
      value: current === null ? EVERYONE : SOME, nullable: false,
      onChanged: (v) => {
        some.value = v === SOME;
        access.setVisibility(table, remote, v === SOME ? access.visibilityOf(table, remote) ?? [] : null);
      }});
    const known = this._known(current ?? []);
    const chips = new ChipsInput({label: 'Groups', name: 'visibleTo', items: known.map(ManifestContextPanel._item),
      value: (current ?? []).map((g) => g.id), enabled: some, emptyText: 'No groups to pick from',
      onChanged: (ids) => access.setVisibility(table, remote, ids.map((id) => known.find((g) => g.id === id)!))});
    const form = new Form({layout: 'wide'});
    form.add(radio).add(chips);
    // a pick becomes a pressed chip; the chips stay the place to unpress it
    if (this._picker !== undefined) {
      const picker = this._picker((principal) => {
        if (!known.some((g) => g.id === principal.id)) {
          known.push(principal);
          chips.setItems(known.map(ManifestContextPanel._item));
        }
        if (!chips.value.peek().includes(principal.id))
          chips.value.value = [...chips.value.peek(), principal.id];
      });
      chips.box.append(picker);
      chips.effect(() => picker.style.display = some.value ? '' : 'none');
    }
    section.add(form, ManifestContextPanel._note(access.editing ?
      'A restricted column is absent from every other user\'s reads, filters and forms; a group let in gets View, ' +
      'and Edit where it may edit the table.' :
      'Applied after Create as a column restriction; a restricted column is absent from every other user\'s ' +
      'reads, filters and forms.'));
    return section;
  }

  private _relation(r: RelationView, caption: string): HTMLElement {
    const row = div([span(caption), badge(r.ref ? 'ref' : 'plain value', {variant: r.ref ? 'accent' : 'warning'}),
      span(r.reason, 'u2-manifest-relation-why')], 'u2-manifest-relation');
    if (!this.offer.editable)
      return row;
    if (r.canFix)
      row.append(link(`include ${r.targetTable}`, () => this.model.includeTable(r.targetTable, true)));
    else if (r.suggested)
      row.append(link('make it a ref', () => this.model.setRef(r.table, r.column, true)));
    else if (r.promoted)
      row.append(link('keep it a plain value', () => this.model.setRef(r.table, r.column, false)));
    return row;
  }

  /** The grid over [rows], handing every change to [onChanged] with the groups resolved;
   * `inheritedGroups` label the inherited rows, which carry group ids alone. */
  private _accessGrid(name: string, rows: Omit<AccessGrant, 'scope'>[], inherited: InheritedAccessRow[],
    locked: string[], onChanged: (rows: Omit<AccessGrant, 'scope'>[]) => void,
    options: {inheritedGroups?: AccessPrincipal[], enabled?: boolean} = {}): AccessGrid {
    const known = this._known([...rows.map((g) => g.group), ...(options.inheritedGroups ?? [])]);
    const picker = this._picker;
    return new AccessGrid({name: `access-${name}`, inline: true,
      capabilities: CAPABILITIES, principals: known.map(ManifestContextPanel._item), inherited, locked,
      enabled: this.offer.editable && options.enabled !== false, value: rows.map((g) => ManifestContextPanel._row(g)),
      picker: picker === undefined ? undefined : (add) => picker((principal) => {
        if (!known.some((g) => g.id === principal.id))
          known.push(principal);
        add({value: principal.id, label: principal.label});
      }),
      onChanged: (changed) => onChanged(changed.map((r) => ({
        group: known.find((g) => g.id === r.principal)!,
        view: r.can.view === true, edit: r.can.edit === true, delete: r.can.delete === true})))});
  }

  /** The offered groups plus the ones already granted — a row reopened from elsewhere keeps its label. */
  private _known(granted: AccessPrincipal[]): AccessPrincipal[] {
    return [...this._groups, ...granted.filter((g) => !this._groups.some((k) => k.id === g.id))];
  }

  private static _item(g: AccessPrincipal): {value: string, label: string} {
    return {value: g.id, label: g.label};
  }

  /** An editor where the offer allows an edit — and the field is not `locked` for this item — the
   * value as text otherwise; the hint sits on the same line, after the editor — the input's
   * postfix, or a span after the text. */
  private _field(form: Form, label: string, name: string, text: string,
    build: (hint: {postfix?: string}) => Input<any>, hint?: string, editable = true): void {
    if (this.offer.editable && editable)
      form.add(build(hint === undefined ? {} : {postfix: hint}));
    else {
      const field = ObjectForm.readonlyField(label, name, text);
      if (hint !== undefined)
        field.value.append(ManifestContextPanel._hint(hint));
      form.addElement(field.row);
    }
  }

  private static _hint(text: string): HTMLElement {
    return span(text, 'u2-manifest-panel-hint');
  }

  private static _row(g: Omit<AccessGrant, 'scope'>): AccessRow {
    return {principal: g.group.id, can: {view: g.view, edit: g.edit, delete: g.delete}};
  }

  /** The catalog's verdict on a registered item, on its title, with its reason. */
  private static _drift(title: HTMLElement, drift: DriftView | undefined): void {
    const label = ManifestTree.driftLabel(drift);
    if (label !== null)
      title.append(badge(label, {variant: drift!.blocks ? 'error' : 'warning'}), span(drift!.reason, 'u2-manifest-node-hint'));
  }

  private static _creator(): InheritedAccessRow {
    return {principal: CREATOR, can: {view: true, edit: true, delete: true}, from: ''};
  }

  private static _title(kind: string, name: string): HTMLElement {
    const el = div([span(kind), span(name, 'u2-manifest-node-hint')], 'u2-manifest-panel-title');
    el.dataset.u2Part = 'title';
    return el;
  }

  private static _note(text: string): HTMLElement {
    return div([span(text)], 'u2-manifest-panel-note');
  }
}

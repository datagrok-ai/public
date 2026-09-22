/* The Design step's right pane: one panel per selection kind, rebuilt under a fresh scope when
   the selection or the model changes (the context-panel pattern of `propertyEditor`). Schema:
   name, friendly name, writable, the access grid for every table. Table: logical and friendly
   names, the key (read-only), the name and searchable pickers, the read-only opt-out under a
   writable binding, the relationships with an "include <target>" link where including the
   target would make a ref, the access grid with the schema's rows inherited. Column: logical
   name, the mapped type beside the warehouse's, required, name, searchable, its reference with
   its status, visibility. The field offer follows the editor context: `view` renders every
   field as text, never as a dead input. Access is modelled beside the manifest and rendered
   here; the dialog applies it after Create. */
import {Control} from '../../../core/component.js';
import {Scope} from '../../../core/scope.js';
import {signal, ReadonlySignal} from '../../../core/signals.js';
import {div, link, span} from '../../../core/elements.js';
import type {Input} from '../../../core/input-base.js';
import {Form} from '../../../components/forms/form.js';
import {Section} from '../../../components/containers/section.js';
import {TextInput} from '../../../components/inputs/text-input.js';
import {BoolInput} from '../../../components/inputs/bool-input.js';
import {ChoiceInput} from '../../../components/inputs/choice-input.js';
import {RadioInput} from '../../../components/inputs/radio-input.js';
import {ChipsInput} from '../../../components/inputs/chips-input.js';
import {AccessGrid} from '../../../components/forms/access-grid.js';
import type {AccessRow, InheritedAccessRow} from '../../../components/forms/access-grid.js';
import {badge} from '../../../components/display/badge.js';
import {ObjectForm} from '../../forms/object-form.js';
import type {EditorContext, FieldOffer} from './editor-context.js';
import {fieldOffer} from './editor-context.js';
import {AccessModel, ManifestModel} from './manifest-model.js';
import type {AccessGrant, AccessPrincipal, AccessScope, ColumnView, ManifestDiagnostic, ManifestSelection,
  RelationView} from './manifest-model.js';
import {ManifestTree} from './manifest-tree.js';

export interface ManifestContextPanelOptions {
  selected: ReadonlySignal<ManifestSelection>;
  access: AccessModel;
  context: EditorContext;
  /** The groups the access pickers offer. */
  groups?: AccessPrincipal[];
  /** Diagnostics addressed by manifest path; the ones about the selected node are listed on top. */
  diagnostics?: ReadonlySignal<ManifestDiagnostic[]>;
  /** Why the Writable switch cannot be turned on here; the switch is then disabled with this hint. */
  writableDisabled?: string;
}

const SCHEMA: AccessScope = {kind: 'schema'};

const CAPABILITIES = [{name: 'view', label: 'View'}, {name: 'edit', label: 'Edit'}, {name: 'delete', label: 'Delete'}];
const CREATOR = 'You (creator)';
const EVERYONE = 'Everyone who may see the row';
const SOME = 'Only these groups';
const NONE = '— none —';

export class ManifestContextPanel extends Control {
  readonly offer: FieldOffer;

  private readonly _access: AccessModel;
  private readonly _groups: AccessPrincipal[];
  private readonly _diagnostics: ReadonlySignal<ManifestDiagnostic[]> | undefined;
  private readonly _writableDisabled: string | undefined;
  private _shown: Scope | undefined;

  constructor(readonly model: ManifestModel, options: ManifestContextPanelOptions) {
    super();
    this.offer = fieldOffer(options.context);
    this._access = options.access;
    this._groups = options.groups ?? [];
    this._diagnostics = options.diagnostics;
    this._writableDisabled = options.writableDisabled;
    this.root.classList.add('u2-manifest-panel');
    this.root.dataset.u2 = 'manifest-panel';
    this.own(() => this._shown?.dispose());
    this.effect(() => {
      const selection = options.selected.value;
      this.model.revision.value;
      this.model.name.value;
      this._shown?.dispose();
      const scope = new Scope();
      this._shown = scope;
      Scope.runWith(scope, () => this._show(selection));
    });
  }

  private _show(selection: ManifestSelection): void {
    const problems = (this._diagnostics?.value ?? [])
      .filter((d) => ManifestModel.sameSelection(this.model.resolvePath(d.path), selection));
    const content = selection.kind === 'schema' ? this._schema() :
      selection.kind === 'table' ? this._table(selection.table) :
        this._column(selection.table, selection.column);
    if (problems.length > 0) {
      content.splice(1, 0, div(problems.map((p) =>
        div([badge(p.code, {variant: 'error'}), span(p.message)], 'u2-manifest-panel-diagnostic')),
      'u2-manifest-panel-diagnostics'));
    }
    this.root.replaceChildren(...content.map((c) => Control.is(c) ? c.root : c));
  }

  private _schema(): (HTMLElement | Control)[] {
    const model = this.model;
    const form = new Form({layout: 'wide'});
    const name = model.name.peek();
    this._field(form, 'Name', 'name', name, () => {
      const input = new TextInput({label: 'Name', name: 'name', value: name, commitOn: 'change',
        onChanged: (v) => model.name.value = v});
      input.addValidator((v) => model.checkSchemaName(v));
      return input;
    });
    form.addElement(ManifestContextPanel._note(
      `registered as ext_${name}; lowercase letters, digits, underscores`));
    const friendly = model.friendlyName.peek();
    this._field(form, 'Friendly name', 'friendlyName', friendly, () => new TextInput({label: 'Friendly name',
      name: 'friendlyName', value: friendly, commitOn: 'change', onChanged: (v) => model.friendlyName.value = v}));
    const writable = model.writable.peek();
    if (this.offer.writable) {
      const hint = this._writableDisabled;
      this._field(form, 'Writable', 'writable', writable ? 'Yes' : 'No', () => new BoolInput({label: 'Writable',
        name: 'writable', value: writable, enabled: hint === undefined, tooltipText: hint,
        onChanged: (v) => model.setWritable(v)}));
      form.addElement(ManifestContextPanel._note(hint === undefined ?
        '— users with Edit on a table may insert, update and delete rows in the warehouse' : `— ${hint}`));
    }
    form.addElement(ManifestContextPanel._note('Queries run as the platform service with the connection\'s ' +
      'stored credentials; users need View on a table, nothing on the connection.'));
    const access = new Section({title: 'Access — every table', collapsible: false});
    access.add(this._accessGrid(SCHEMA, [ManifestContextPanel._creator()], writable ? [] : ['edit', 'delete']));
    access.add(ManifestContextPanel._note('Applied after Create as the same grant on every included table — ' +
      'there is no schema-wide row access. Edit and Delete need a writable binding.'));
    return [ManifestContextPanel._title('Schema', name), form, access];
  }

  private _table(remote: string): (HTMLElement | Control)[] {
    const model = this.model;
    const table = model.table(remote);
    if (table === undefined)
      return [ManifestContextPanel._title('Table', remote)];
    const title = ManifestContextPanel._title('Table', remote);
    if (!table.included) {
      const status = table.bindable ? 'not included' : ManifestTree.shortReason(table.code);
      title.append(badge(status, {variant: table.bindable ? 'warning' : 'error'}),
        span(table.reason ?? 'check it in the tree to expose it', 'u2-manifest-node-hint'));
      if (!table.bindable)
        return [title];
    }
    const columns = model.columns(remote).peek();
    const form = new Form({layout: 'wide'});
    this._field(form, 'Logical name', 'logical', table.logical, () => {
      const input = new TextInput({label: 'Logical name', name: 'logical', value: table.logical, commitOn: 'change',
        onChanged: (v) => model.renameTable(remote, v)});
      input.addValidator((v) => model.checkTableName(remote, v));
      return input;
    });
    this._field(form, 'Friendly name', 'friendlyName', table.friendlyName, () => new TextInput({label: 'Friendly name',
      name: 'friendlyName', value: table.friendlyName, commitOn: 'change',
      onChanged: (v) => model.setFriendlyName(remote, v)}));
    const key = table.key.map((k) => columns.find((c) => c.remote === k)?.logical ?? k).join(', ');
    form.addElement(ObjectForm.readonlyField('Key', 'key', `${key} — the primary key; the row id encodes it`).row);
    const strings = columns.filter((c) => c.supported && c.included && c.type === 'string')
      .map((c) => ({value: c.remote, label: c.logical}));
    const pick = (label: string, name: string, current: ColumnView | undefined,
      apply: (remote: string | null) => void): void => {
      this._field(form, label, name, current?.logical ?? NONE, () => new ChoiceInput({label, name, items: strings,
        value: current?.remote ?? null, onChanged: apply}));
    };
    if (this.offer.nameColumn)
      pick('Name column', 'nameColumn', columns.find((c) => c.isName), (r) => model.setNameColumn(remote, r));
    if (this.offer.searchable)
      pick('Searchable', 'searchable', columns.find((c) => c.searchable), (r) => model.setSearchable(remote, r));
    if (this.offer.writable && model.writable.peek()) {
      this._field(form, 'Read-only', 'readOnly', table.readOnly ? 'Yes' : 'No', () => new BoolInput({label: 'Read-only',
        name: 'readOnly', value: table.readOnly, onChanged: (v) => model.setReadOnly(remote, v)}));
      form.addElement(ManifestContextPanel._note('— opt this table out of writes'));
    }
    const relations = new Section({title: 'Relationships', collapsible: false});
    const rels = model.relations.peek().filter((r) => r.table === remote);
    if (rels.length === 0)
      relations.add(ManifestContextPanel._note('The warehouse reports no foreign keys on this table.'));
    for (const r of rels)
      relations.add(this._relation(r, `${r.column} → ${r.targetTable}`));
    const schemaRows = this._access.grantsOf(SCHEMA).peek()
      .map((g): InheritedAccessRow => ({...ManifestContextPanel._row(g), from: '(every table)'}));
    const locked = !model.writable.peek() || table.readOnly ? ['edit', 'delete'] : [];
    const access = new Section({title: 'Access — this table', collapsible: false});
    access.add(this._accessGrid({kind: 'table', table: remote}, [ManifestContextPanel._creator(), ...schemaRows],
      locked));
    return [title, form, relations, access];
  }

  private _column(table: string, remote: string): (HTMLElement | Control)[] {
    const model = this.model;
    const column = model.column(table, remote);
    const title = ManifestContextPanel._title('Column', `${table}.${remote}`);
    if (column === undefined)
      return [title];
    if (!column.supported) {
      title.append(badge('not bindable', {variant: 'error'}), span(column.reason ?? '', 'u2-manifest-node-hint'));
      return [title];
    }
    const form = new Form({layout: 'wide'});
    this._field(form, 'Logical name', 'logical', column.logical, () => {
      const input = new TextInput({label: 'Logical name', name: 'logical', value: column.logical, commitOn: 'change',
        onChanged: (v) => model.renameColumn(table, remote, v)});
      input.addValidator((v) => model.checkColumnName(table, remote, v));
      return input;
    });
    const type = column.type === 'ref' ? `ref → ${column.relation?.targetLogical ?? ''}` : column.type;
    form.addElement(ObjectForm.readonlyField('Type', 'type',
      column.dbType === undefined ? type : `${type} ← ${column.dbType}`).row);
    if (this.offer.required) {
      this._field(form, 'Required', 'required', column.required ? 'Yes' : 'No', () => new BoolInput({label: 'Required',
        name: 'required', value: column.required, enabled: !column.isKey,
        onChanged: (v) => model.setRequired(table, remote, v)}));
      form.addElement(ManifestContextPanel._note(column.isKey ? '— key columns are always required' :
        '— NOT NULL is not reported by the warehouse; tick it if you know'));
    }
    const isString = column.type === 'string';
    if (this.offer.nameColumn) {
      this._field(form, 'Name column', 'isName', column.isName ? 'Yes' : 'No', () => new BoolInput({
        label: 'Name column', name: 'isName', value: column.isName, enabled: isString,
        onChanged: (v) => model.setNameColumn(table, v ? remote : null)}));
      form.addElement(ManifestContextPanel._note('— shown wherever a row is referred to'));
    }
    if (this.offer.searchable) {
      this._field(form, 'Searchable', 'searchable', column.searchable ? 'Yes' : 'No', () => new BoolInput({
        label: 'Searchable', name: 'searchable', value: column.searchable, enabled: isString,
        onChanged: (v) => model.setSearchable(table, v ? remote : null)}));
      form.addElement(ManifestContextPanel._note('— one per table'));
    }
    const content: (HTMLElement | Control)[] = [title, form];
    const relation = column.relation;
    if (relation !== undefined) {
      const reference = new Section({title: 'Reference', collapsible: false});
      reference.add(this._relation(relation, `foreign key to ${relation.targetTable}`));
      content.push(reference);
    }
    content.push(this._visibility(table, remote, column));
    return content;
  }

  private _visibility(table: string, remote: string, column: ColumnView): Section {
    const section = new Section({title: 'Visibility', collapsible: false});
    if (column.isKey) {
      section.add(ManifestContextPanel._note(
        'A key column is visible to everyone who may see the row: the row id encodes it.'));
      return section;
    }
    const access = this._access;
    const current = access.visibilityOf(table, remote);
    if (!this.offer.editable) {
      section.add(ManifestContextPanel._note(current === null ? EVERYONE :
        `${SOME}: ${current.length === 0 ? 'none yet' : current.map((g) => g.label).join(', ')}`));
      return section;
    }
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
    section.add(form, ManifestContextPanel._note('Applied after Create as a column restriction; a restricted ' +
      'column is absent from every other user\'s reads, filters and forms.'));
    return section;
  }

  private _relation(r: RelationView, caption: string): HTMLElement {
    const row = div([span(caption), badge(r.ref ? 'ref' : 'plain value', {variant: r.ref ? 'accent' : 'warning'}),
      span(r.reason, 'u2-manifest-relation-why')], 'u2-manifest-relation');
    if (r.canFix && this.offer.editable)
      row.append(link(`include ${r.targetTable}`, () => this.model.includeTable(r.targetTable, true)));
    return row;
  }

  private _accessGrid(scope: AccessScope, inherited: InheritedAccessRow[], locked: string[]): AccessGrid {
    const access = this._access;
    const rows = access.grantsOf(scope).peek();
    const known = this._known(rows.map((g) => g.group));
    const grid = new AccessGrid({name: `access-${scope.kind === 'schema' ? 'schema' : scope.table}`, inline: true,
      capabilities: CAPABILITIES, principals: known.map(ManifestContextPanel._item), inherited, locked,
      enabled: this.offer.editable, value: rows.map((g) => ManifestContextPanel._row(g)),
      onChanged: (changed) => access.setGrants(scope, changed.map((r) => ({
        group: known.find((g) => g.id === r.principal)!,
        view: r.can.view === true, edit: r.can.edit === true, delete: r.can.delete === true})))});
    return grid;
  }

  /** The offered groups plus the ones already granted — a row reopened from elsewhere keeps its label. */
  private _known(granted: AccessPrincipal[]): AccessPrincipal[] {
    return [...this._groups, ...granted.filter((g) => !this._groups.some((k) => k.id === g.id))];
  }

  private static _item(g: AccessPrincipal): {value: string, label: string} {
    return {value: g.id, label: g.label};
  }

  /** An editor where the offer allows an edit, the value as text otherwise. */
  private _field(form: Form, label: string, name: string, text: string, build: () => Input<any>): void {
    if (this.offer.editable)
      form.add(build());
    else
      form.addElement(ObjectForm.readonlyField(label, name, text).row);
  }

  private static _row(g: AccessGrant): AccessRow {
    return {principal: g.group.id, can: {view: g.view, edit: g.edit, delete: g.delete}};
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

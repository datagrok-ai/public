import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {computed, div, divV, span, plural, statusLine, signal, Signal, Dialog, BoolInput, notify}
  from '@datagrok-libraries/u2';
import {errorCode} from '../../core/domain-backend.js';
import {DomainErrors} from '../errors.js';
import {BindingWizard} from './binding-dialog.js';
import {ManifestEditor} from './manifest-editor.js';
import type {ManifestChange, ManifestJson, ManifestOwnership} from './manifest-model.js';

const OPTIONS = ['friendlyName', 'description', 'choices', 'min', 'max', 'format'] as const;

interface Consequence {
  verb: string;
  noun: [one: string, many: string];
  text: string;
  said: boolean;
  count?: number;
  notes?: HTMLElement[];
}

export class PropertyDialog {
  private static readonly _extensible = new Map<string, Promise<Set<string>>>();

  static readonly UNSUPPORTED = 'Edit properties is unavailable on this server version: it cannot report which ' +
    'columns you added, and saving without that could delete them';

  static extensible(schema: string): Promise<Set<string>> {
    let tables = PropertyDialog._extensible.get(schema);
    if (tables === undefined)
      PropertyDialog._extensible.set(schema, tables = PropertyDialog._readExtensible(schema).catch(() => new Set<string>()));
    return tables;
  }

  static async open(schema: string, table?: string): Promise<boolean> {
    if ((await PropertyDialog.extensible(schema)).size === 0)
      throw new Error(`You may not add properties to ${schema}`);
    const handle = grok.dapi.domains.schema(schema);
    const {manifest, ownership} = await PropertyDialog._baseline(handle);
    const name = table === undefined ? manifest.name : manifest.tables[table]?.['singularName'] as string ?? table;
    const dialog = Dialog.create(`${name.charAt(0).toUpperCase()}${name.slice(1)} properties`, {name: 'properties'});
    const editor = dialog.build(() => new ManifestEditor(null, {context: {mode: 'extend', storage: 'domain'},
      baseline: manifest, ownership, columnOptions: OPTIONS.map((option) => DG.Property.propertyOptions[option]!)}));
    const changes = computed(() => {
      editor.model.revision.value;
      return editor.model.changes();
    });
    const reloaded = signal<string | null>(null);
    let saved = false;
    return new Promise((resolve) => {
      dialog.add(() => divV([statusLine(reloaded), editor.root], 'u2-manifest-properties'))
        .okDisabled(computed(() => changes.value.length === 0 ? 'No changes to save' : null))
        .onOK(async () => saved = await PropertyDialog._save(handle, editor, reloaded), {text: 'SAVE', primary: true})
        .onClose(() => resolve(saved))
        .show({width: 980, height: 640});
      if (table !== undefined)
        void editor.tree.reveal(table);
    });
  }

  private static async _readExtensible(schema: string): Promise<Set<string>> {
    const entity = await PropertyDialog._schema(grok.dapi.domains.schemas, schema);
    if (entity?.managedBy !== 'package' || !await grok.dapi.permissions.check(entity, 'Extend'))
      return new Set();
    const tables = (await grok.dapi.domains.schema(schema).manifest()).tables;
    return new Set(Object.keys(tables).filter((t) => tables[t].extensible === true));
  }

  private static async _baseline(handle: DG.DomainSchemaClient): Promise<{manifest: ManifestJson,
    ownership: ManifestOwnership}> {
    const {extVersion} = await handle.apply({}, {dryRun: true});
    if (extVersion === undefined)
      throw new Error(`${handle.name} answered no extension version: properties cannot be saved safely`);
    const manifest = await (handle.manifest() as Promise<ManifestJson>);
    const ownership = await PropertyDialog._ownership(handle.name, Object.keys(manifest.tables));
    return {manifest: {...manifest, extVersion: String(extVersion)}, ownership};
  }

  private static async _schema(schemas: DG.HttpDataSource<DG.DomainSchema>, name: string):
    Promise<DG.DomainSchema | undefined> {
    return (await schemas.filter(`name = ${JSON.stringify(name)}`).list()).find((s) => s.name === name);
  }

  private static async _ownership(schema: string, tables: string[]): Promise<ManifestOwnership> {
    const entity = await PropertyDialog._schema(grok.dapi.domains.schemas.include('tables.columns.property'), schema);
    if (entity === undefined)
      throw new Error(`${schema} is not a registered schema`);
    const ownership: ManifestOwnership = {};
    for (const t of entity.tables) {
      const {origin, userColumns} = t;
      if ((origin !== 'package' && origin !== 'user') || !Array.isArray(userColumns) ||
        !userColumns.every((c) => typeof c === 'string'))
        throw new Error(PropertyDialog.UNSUPPORTED);
      ownership[t.name] = {origin, userColumns};
    }
    const missing = tables.find((t) => ownership[t] === undefined);
    if (missing !== undefined)
      throw new Error(`${schema}.${missing} is not in the registry yet: reload and try again`);
    return ownership;
  }

  private static async _save(handle: DG.DomainSchemaClient, editor: ManifestEditor,
    reloaded: Signal<string | null>): Promise<boolean> {
    const {payload, changes} = editor.extendPlan();
    try {
      const plan = await handle.apply(payload, {dryRun: true}).catch((e) => {
        const demanded = (e as {body?: {plan?: DG.DomainApplyPlan}}).body?.plan;
        if (demanded === undefined)
          throw e;
        return demanded;
      });
      if (plan.destructive && !await PropertyDialog._confirm(changes, plan))
        return false;
      await handle.apply({...payload, confirmDestructive: plan.destructive});
    }
    catch (e) {
      if (errorCode(e) !== 'version-conflict') {
        editor.diagnostics.value = BindingWizard.issues(e);
        throw new Error(DomainErrors.message(e));
      }
      const {manifest, ownership} = await PropertyDialog._baseline(handle);
      const report = editor.rebase(manifest, null, null, {friendlyName: '', description: '', ownership});
      const lost = [...report.conflicts, ...report.dropped].map((op) => op.note ?? `${op.text} no longer applies`);
      reloaded.value = 'Someone else changed these properties meanwhile; reloaded with your edits on top. ' +
        `${lost.length > 0 ? `${lost.join('; ')}. ` : ''}Review and save again.`;
      return false;
    }
    grok.dapi.domains.invalidateUiCaches();
    notify.info(`Saved ${plural(changes.length, 'change', 'changes')}`);
    return true;
  }

  private static _confirm(changes: ManifestChange[], plan: DG.DomainApplyPlan): Promise<boolean> {
    const counted = (n: number | undefined, one: string, many: string): HTMLElement[] =>
      n !== undefined && n > 0 ? [span(`${plural(n, one, many)} deleted`, 'u2-binding-change-note')] : [];
    const removed = (table: string, column?: string): boolean =>
      changes.some((c) => c.removes && c.table === table && c.column === column);
    const edited = (table: string, column: string): boolean => changes.some((c) => c.table === table &&
      c.column === column && (c.id.endsWith(':decl') || c.id.endsWith(':ref')));
    const {drops: {tables, columns}, typeChanges, uniqueDrops, autoNumberChanges, businessKeyChanges} = plan;
    const constraints = Object.entries(plan.constraintDrops);
    const values = (d: {nonNullValues?: number}): number => d.nonNullValues ?? 0;
    const rows = (t: {liveRows?: number}): number => t.liveRows ?? 0;
    const column = (t: string, c: string): string => `Column ${t}.${c}`;
    const consequences: Consequence[] = [
      ...columns.filter((d) => values(d) > 0).map((d): Consequence => ({verb: 'delete', noun: ['stored value',
        'stored values'], count: values(d), text: `${column(d.table, d.column)} dropped`,
      said: removed(d.table, d.column), notes: counted(values(d), 'value', 'values')})),
      ...tables.filter((t) => rows(t) > 0).map((t): Consequence => ({verb: 'delete', noun: ['row', 'rows'],
        count: rows(t), text: `Table ${t.name} dropped`, said: removed(t.name),
        notes: counted(rows(t), 'row', 'rows')})),
      ...columns.filter((d) => values(d) === 0).map((d): Consequence => ({verb: 'drop', noun: ['column', 'columns'],
        text: `${column(d.table, d.column)} dropped`, said: removed(d.table, d.column)})),
      ...tables.filter((t) => rows(t) === 0).map((t): Consequence => ({verb: 'drop', noun: ['table', 'tables'],
        text: `Table ${t.name} dropped`, said: removed(t.name)})),
      ...typeChanges.map((t): Consequence => ({verb: 'change the type of', noun: ['column', 'columns'],
        text: `${column(t.table, t.column)}: type changes to ${t.to}`, said: edited(t.table, t.column)})),
      ...uniqueDrops.map((u): Consequence => ({verb: 'drop', noun: ['unique constraint', 'unique constraints'],
        text: `${column(u.table, u.column)}: unique constraint dropped`, said: edited(u.table, u.column)})),
      ...Object.entries(autoNumberChanges).flatMap(([t, cs]) => cs.map((c): Consequence => ({
        verb: 'change the auto-numbering of', noun: ['column', 'columns'],
        text: `${column(t, c)}: auto-numbering changes`, said: edited(t, c)}))),
      ...businessKeyChanges.map((k): Consequence => ({verb: 'change', noun: ['business key', 'business keys'],
        text: `Table ${k.table}: business key ${k.businessKey.length === 0 ? 'removed' :
          `changes to ${k.businessKey.join(', ')}`}`, said: false})),
      ...constraints.flatMap(([table, names]) => names.map((k): Consequence => ({verb: 'drop',
        noun: ['constraint', 'constraints'], text: `Table ${table}: constraint ${k} dropped`, said: false}))),
    ];
    const lost = (c: ManifestChange): HTMLElement[] => !c.removes ? [] : c.column === undefined ?
      counted(tables.find((t) => t.name === c.table)?.liveRows, 'row', 'rows') :
      counted(columns.find((d) => d.table === c.table && d.column === c.column)?.nonNullValues, 'value', 'values');
    const lines = [
      ...changes.map((c) => div([span(c.text), ...lost(c)],
        `u2-binding-change${c.removes ? ' u2-binding-change-removes' : ''}`)),
      ...consequences.filter((c) => !c.said).map((c) =>
        div([span(c.text), ...c.notes ?? []], 'u2-binding-change u2-binding-change-removes')),
    ];
    const totals = new Map<string, {verb: string, noun: [string, string], count: number}>();
    for (const {verb, noun, count = 1} of consequences) {
      const key = `${verb} ${noun[1]}`;
      totals.set(key, {verb, noun, count: (totals.get(key)?.count ?? 0) + count});
    }
    const parts = [...totals.values()].map(({verb, noun: [one, many], count}) => `${verb} ${plural(count, one, many)}`);
    const words = parts.length === 0 ? 'apply the destructive changes' : parts.length === 1 ? parts[0] :
      `${parts.slice(0, -1).join(', ')} and ${parts[parts.length - 1]}`;
    const confirmed = signal(false);
    return new Promise((resolve) => {
      let answer = false;
      const dialog = Dialog.create('Save these changes?', {name: 'properties-confirm'});
      const box = dialog.build(() => new BoolInput({label: words.charAt(0).toUpperCase() + words.slice(1),
        name: 'confirmDestructive', onChanged: (v) => confirmed.value = v}));
      dialog.add(divV([...lines, box.root], 'u2-binding-changes'))
        .okDisabled(computed(() => confirmed.value ? null : 'Confirm the destructive changes first'))
        .onOK(() => answer = true, {text: 'SAVE', primary: true})
        .onClose(() => resolve(answer))
        .show({width: 520});
    });
  }
}

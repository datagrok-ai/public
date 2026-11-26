import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {computed, div, divV, span, plural, statusLine, signal, Signal, Dialog, BoolInput, notify}
  from '@datagrok-libraries/u2';
import {DomainErrors} from '../errors.js';
import {BindingWizard} from './binding-dialog.js';
import {ManifestEditor} from './manifest-editor.js';
import type {ManifestChange, ManifestJson, ManifestOwnership} from './manifest-model.js';

const OPTIONS = ['friendlyName', 'description', 'choices', 'min', 'max', 'format'] as const;

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

  // the token is read FIRST, so a change landing between the reads is a conflict, never an overwrite
  private static async _baseline(handle: DG.DomainSchemaClient): Promise<{manifest: ManifestJson,
    ownership: ManifestOwnership}> {
    const {extVersion} = await handle.apply({}, {dryRun: true});
    if (extVersion === undefined)
      throw new Error(`${handle.name} answered no extension version: properties cannot be saved safely`);
    const manifest = await (handle.manifest() as Promise<ManifestJson>);
    const ownership = await PropertyDialog._ownership(handle.name, Object.keys(manifest.tables));
    return {manifest: {...manifest, extVersion: String(extVersion)}, ownership};
  }

  // the filter may match names in another case: the schema is the one spelled exactly so
  private static async _schema(schemas: DG.HttpDataSource<DG.DomainSchema>, name: string):
    Promise<DG.DomainSchema | undefined> {
    return (await schemas.filter(`name = ${JSON.stringify(name)}`).list()).find((s) => s.name === name);
  }

  // a column read as the plugin's is left out of the full-state extension, which the server takes as a drop
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
      if (DomainErrors.codeOf(e) !== 'version-conflict') {
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

  // built from the server's plan, not the client's diff: SAVE confirms the whole plan, each part once
  private static _confirm(changes: ManifestChange[], plan: DG.DomainApplyPlan): Promise<boolean> {
    const counted = (n: number | undefined, one: string, many: string): HTMLElement[] =>
      n !== undefined && n > 0 ? [span(`${plural(n, one, many)} deleted`, 'u2-binding-change-note')] : [];
    const removed = (table: string, column?: string): boolean =>
      changes.some((c) => c.removes && c.table === table && c.column === column);
    const edited = (table: string, column: string): boolean => changes.some((c) => c.table === table &&
      c.column === column && (c.id.endsWith(':decl') || c.id.endsWith(':ref')));
    const {tables, columns} = plan.drops;
    const constraints = Object.entries(plan.constraintDrops as {[table: string]: string[]} | null ?? {});
    const values = (d: DG.DomainApplyPlan['drops']['columns'][number]) => d.nonNullValues ?? 0;
    const rows = (t: DG.DomainApplyPlan['drops']['tables'][number]) => t.liveRows ?? 0;
    const column = (t: string, c: string) => `Column ${t}.${c}`;
    type Item = {n: number, text: string, said: boolean, notes?: HTMLElement[]};
    const item = (text: string, said: boolean, n = 1, notes?: HTMLElement[]): Item => ({n, text, said, notes});
    // per kind: what the plan does, a line each unless a client change says it already, and the consent's words
    const kinds: [verb: string, one: string, many: string, items: Item[]][] = [
      ['delete', 'stored value', 'stored values', columns.filter((d) => values(d) > 0).map((d) =>
        item(`${column(d.table, d.column)} dropped`, removed(d.table, d.column), values(d),
          counted(values(d), 'value', 'values')))],
      ['delete', 'row', 'rows', tables.filter((t) => rows(t) > 0).map((t) =>
        item(`Table ${t.name} dropped`, removed(t.name), rows(t), counted(rows(t), 'row', 'rows')))],
      ['drop', 'column', 'columns', columns.filter((d) => values(d) === 0).map((d) =>
        item(`${column(d.table, d.column)} dropped`, removed(d.table, d.column)))],
      ['drop', 'table', 'tables', tables.filter((t) => rows(t) === 0).map((t) =>
        item(`Table ${t.name} dropped`, removed(t.name)))],
      ['change the type of', 'column', 'columns', (plan.typeChanges ?? []).map((t) =>
        item(`${column(t.table, t.column)}: type changes to ${t.to}`, edited(t.table, t.column)))],
      ['drop', 'unique constraint', 'unique constraints', (plan.uniqueDrops ?? []).map((u) =>
        item(`${column(u.table, u.column)}: unique constraint dropped`, edited(u.table, u.column)))],
      ['change the auto-numbering of', 'column', 'columns', Object.entries(plan.autoNumberChanges ?? {})
        .flatMap(([t, cs]) => cs.map((c) => item(`${column(t, c)}: auto-numbering changes`, edited(t, c))))],
      ['change', 'business key', 'business keys', (plan.businessKeyChanges ?? []).map((k) => item(`Table ${k.table}: ` +
        `business key ${k.businessKey.length === 0 ? 'removed' : `changes to ${k.businessKey.join(', ')}`}`, false))],
      ['drop', 'constraint', 'constraints', constraints.flatMap(([table, names]) =>
        names.map((k) => item(`Table ${table}: constraint ${k} dropped`, false)))],
    ];
    const lost = (c: ManifestChange): HTMLElement[] => !c.removes ? [] : c.column === undefined ?
      counted(tables.find((t) => t.name === c.table)?.liveRows, 'row', 'rows') :
      counted(columns.find((d) => d.table === c.table && d.column === c.column)?.nonNullValues, 'value', 'values');
    const lines = [
      ...changes.map((c) => div([span(c.text), ...lost(c)],
        `u2-binding-change${c.removes ? ' u2-binding-change-removes' : ''}`)),
      ...kinds.flatMap(([, , , items]) => items.filter((i) => !i.said).map((i) =>
        div([span(i.text), ...i.notes ?? []], 'u2-binding-change u2-binding-change-removes'))),
    ];
    const parts = kinds.flatMap(([verb, one, many, items]) => {
      const n = items.reduce((sum, i) => sum + i.n, 0);
      return n === 0 ? [] : [`${verb} ${plural(n, one, many)}`];
    });
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

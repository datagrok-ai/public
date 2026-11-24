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
    const entity = (await grok.dapi.domains.schemas.list()).find((s) => s.name === schema);
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
    const manifest = await handle.manifest() as unknown as ManifestJson;
    const ownership = await PropertyDialog._ownership(handle.name);
    for (const table of Object.keys(manifest.tables)) {
      if (ownership[table] === undefined)
        throw new Error(`${handle.name}.${table} is not in the registry yet: reload and try again`);
    }
    return {manifest: {...manifest, extVersion: String(extVersion)}, ownership};
  }

  // a column read as the plugin's is left out of the full-state extension, which the server takes as a drop
  private static async _ownership(schema: string): Promise<ManifestOwnership> {
    const entity = (await grok.dapi.domains.schemas.include('tables.columns.property').list())
      .find((s) => s.name === schema);
    if (entity === undefined)
      throw new Error(`${schema} is not a registered schema`);
    const ownership: ManifestOwnership = {};
    for (const t of entity.tables) {
      let origin: unknown;
      let userColumns: unknown;
      try {
        origin = t.origin;
        userColumns = t.userColumns;
      } catch {
      }
      if ((origin !== 'package' && origin !== 'user') || !Array.isArray(userColumns) ||
        !userColumns.every((c) => typeof c === 'string'))
        throw new Error(PropertyDialog.UNSUPPORTED);
      ownership[t.name] = {origin, userColumns: userColumns as string[]};
    }
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

  // built from the server's plan, not the client's diff: SAVE confirms the whole plan
  private static _confirm(changes: ManifestChange[], plan: DG.DomainApplyPlan): Promise<boolean> {
    const counted = (n: number | undefined, one: string, many: string): HTMLElement[] =>
      n !== undefined && n > 0 ? [span(`${plural(n, one, many)} deleted`, 'u2-binding-change-note')] : [];
    const removed = (table: string, column?: string): ManifestChange | undefined =>
      changes.find((c) => c.removes && c.table === table && c.column === column);
    const edited = (table: string, column: string): boolean => changes.some((c) => c.table === table &&
      c.column === column && (c.id.endsWith(':decl') || c.id.endsWith(':ref')));
    const line = (text: string, notes: HTMLElement[] = []): HTMLElement =>
      div([span(text), ...notes], 'u2-binding-change u2-binding-change-removes');
    const {tables, columns} = plan.drops;
    const constraints = Object.entries(plan.constraintDrops as {[table: string]: string[]} | null ?? {});
    const autoNumbered = Object.entries(plan.autoNumberChanges ?? {})
      .flatMap(([table, cs]) => cs.map((column) => ({table, column})));
    const lost = (c: ManifestChange): HTMLElement[] => !c.removes ? [] : c.column === undefined ?
      counted(tables.find((t) => t.name === c.table)?.liveRows, 'row', 'rows') :
      counted(columns.find((d) => d.table === c.table && d.column === c.column)?.nonNullValues, 'value', 'values');
    const lines = [
      ...changes.map((c) => div([span(c.text), ...lost(c)],
        `u2-binding-change${c.removes ? ' u2-binding-change-removes' : ''}`)),
      ...tables.filter((t) => removed(t.name) === undefined)
        .map((t) => line(`Table ${t.name} dropped`, counted(t.liveRows, 'row', 'rows'))),
      ...columns.filter((d) => removed(d.table, d.column) === undefined)
        .map((d) => line(`Column ${d.table}.${d.column} dropped`, counted(d.nonNullValues, 'value', 'values'))),
      ...(plan.typeChanges ?? []).filter((t) => !edited(t.table, t.column))
        .map((t) => line(`Column ${t.table}.${t.column}: type changes to ${t.to}`)),
      ...(plan.uniqueDrops ?? []).filter((u) => !edited(u.table, u.column))
        .map((u) => line(`Column ${u.table}.${u.column}: unique constraint dropped`)),
      ...autoNumbered.filter((a) => !edited(a.table, a.column))
        .map((a) => line(`Column ${a.table}.${a.column}: auto-numbering changes`)),
      ...(plan.businessKeyChanges ?? []).map((k) => line(`Table ${k.table}: business key ` +
        (k.businessKey.length === 0 ? 'removed' : `changes to ${k.businessKey.join(', ')}`))),
      ...constraints.flatMap(([table, names]) => names.map((k) => line(`Table ${table}: constraint ${k} dropped`))),
    ];
    const label = PropertyDialog._consent(plan, autoNumbered.length,
      constraints.reduce((n, [, names]) => n + names.length, 0));
    const confirmed = signal(false);
    return new Promise((resolve) => {
      let answer = false;
      const dialog = Dialog.create('Save these changes?', {name: 'properties-confirm'});
      const box = dialog.build(() => new BoolInput({label, name: 'confirmDestructive',
        onChanged: (v) => confirmed.value = v}));
      dialog.add(divV([...lines, box.root], 'u2-binding-changes'))
        .okDisabled(computed(() => confirmed.value ? null : 'Confirm the destructive changes first'))
        .onOK(() => answer = true, {text: 'SAVE', primary: true})
        .onClose(() => resolve(answer))
        .show({width: 520});
    });
  }

  private static _consent(plan: DG.DomainApplyPlan, autoNumbered: number, constraints: number): string {
    const values = plan.drops.columns.reduce((n, d) => n + (d.nonNullValues ?? 0), 0);
    const rows = plan.drops.tables.reduce((n, t) => n + (t.liveRows ?? 0), 0);
    const parts = [
      ...(values > 0 ? [`delete ${plural(values, 'stored value', 'stored values')}`] : []),
      ...(rows > 0 ? [`delete ${plural(rows, 'row', 'rows')}`] : []),
      ...([
        {n: plan.drops.columns.filter((d) => !d.nonNullValues).length, verb: 'drop', one: 'column', many: 'columns'},
        {n: plan.drops.tables.filter((t) => !t.liveRows).length, verb: 'drop', one: 'table', many: 'tables'},
        {n: plan.typeChanges?.length ?? 0, verb: 'change the type of', one: 'column', many: 'columns'},
        {n: plan.uniqueDrops?.length ?? 0, verb: 'drop', one: 'unique constraint', many: 'unique constraints'},
        {n: autoNumbered, verb: 'change the auto-numbering of', one: 'column', many: 'columns'},
        {n: plan.businessKeyChanges?.length ?? 0, verb: 'change', one: 'business key', many: 'business keys'},
        {n: constraints, verb: 'drop', one: 'constraint', many: 'constraints'},
      ]).filter((p) => p.n > 0).map((p) => `${p.verb} ${plural(p.n, p.one, p.many)}`),
    ];
    if (parts.length === 0)
      return 'Apply the destructive changes';
    const text = parts.length === 1 ? parts[0] : `${parts.slice(0, -1).join(', ')} and ${parts[parts.length - 1]}`;
    return text.charAt(0).toUpperCase() + text.slice(1);
  }
}

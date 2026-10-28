import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {computed, div, divV, span, plural, statusLine, signal, Signal, Dialog, BoolInput, notify}
  from '@datagrok-libraries/u2';
import {DomainErrors} from '../errors.js';
import {BindingWizard} from './binding-dialog.js';
import {ManifestEditor} from './manifest-editor.js';
import type {ManifestChange, ManifestJson} from './manifest-model.js';

const OPTIONS = ['friendlyName', 'description', 'choices', 'min', 'max', 'format'] as const;

export class PropertyDialog {
  private static readonly _extensible = new Map<string, Promise<Set<string>>>();

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
    const manifest = await handle.manifest() as unknown as ManifestJson;
    const name = table === undefined ? manifest.name : manifest.tables[table]?.['singularName'] as string ?? table;
    const dialog = Dialog.create(`${name.charAt(0).toUpperCase()}${name.slice(1)} properties`, {name: 'properties'});
    const editor = dialog.build(() => new ManifestEditor(null, {context: {mode: 'extend', storage: 'domain'},
      baseline: manifest, columnOptions: OPTIONS.map((option) => DG.Property.propertyOptions[option]!)}));
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
      const report = editor.rebase(await handle.manifest() as unknown as ManifestJson, null, null,
        {friendlyName: '', description: ''});
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
    const values = (c: ManifestChange): number => !c.removes ? 0 :
      plan.drops.columns.find((d) => d.table === c.table && d.column === c.column)?.nonNullValues ?? 0;
    const lost = changes.reduce((n, c) => n + values(c), 0);
    const confirmed = signal(lost === 0);
    const lines = changes.map((c) => div([span(c.text),
      ...(values(c) > 0 ? [span(`${plural(values(c), 'value', 'values')} deleted`, 'u2-binding-change-note')] : [])],
    `u2-binding-change${c.removes ? ' u2-binding-change-removes' : ''}`));
    return new Promise((resolve) => {
      let answer = false;
      const dialog = Dialog.create('Save these changes?', {name: 'properties-confirm'});
      const box = dialog.build(() => new BoolInput({label: `Delete ${plural(lost, 'stored value', 'stored values')}`,
        name: 'confirmDestructive', visible: lost > 0, onChanged: (v) => confirmed.value = v}));
      dialog.add(divV([...lines, box.root], 'u2-binding-changes'))
        .okDisabled(computed(() => confirmed.value ? null : 'Confirm the deletion first'))
        .onOK(() => answer = true, {text: 'SAVE', primary: true})
        .onClose(() => resolve(answer))
        .show({width: 520});
    });
  }
}

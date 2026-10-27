/* `domains.editProperties` — the user's own columns on a package schema's extensible tables: the
   `ManifestEditor` in extend mode over the registered manifest. SAVE dry-runs the extend section;
   a plan that drops values asks first, naming how many go; a version conflict reloads the schema
   and replays the edits, saying what no longer applies; a refusal lands on the rows it names. The
   open views of a changed table are rebuilt over the new columns. The server owns every rule. */
import * as grok from 'datagrok-api/grok';
import type * as DG from 'datagrok-api/dg';
import {computed, signal, div, divV, span, plural, Dialog, BoolInput, badge, notify} from '@datagrok-libraries/u2';
import {DomainErrors} from '../errors.js';
import {BindingWizard} from './binding-dialog.js';
import {ManifestEditor} from './manifest-editor.js';
import type {ExtendPayload} from './manifest-editor.js';
import type {ManifestChange, ManifestJson, RebaseReport} from './manifest-model.js';

export interface EditPropertiesOptions {
  /** The table to open at — the one a "properties" link stands for. */
  table?: string;
  /** The dialog's title; the table's, else the schema's, "properties". */
  title?: string;
  /** Runs once a save landed, with the tables whose columns changed — what rebuilds their views. */
  onSaved?: (tables: string[]) => void | Promise<void>;
}

/** Opens the editor; resolves to whether a change was saved. Refuses, by name, a schema that is not
 * a package's or a caller without Extend on it. */
export async function editProperties(schema: string, options: EditPropertiesOptions = {}): Promise<boolean> {
  if (!await PropertyDialog.canExtend(schema))
    throw new Error(`You don't have permission to add properties to ${schema}`);
  const handle = grok.dapi.domains.schema(schema);
  const manifest = await handle.manifest();
  if (manifest.extVersion === undefined)
    throw new Error(`${schema} is not a plugin's schema — its tables are edited, not extended`);
  return new PropertyDialog(schema, handle, manifest as unknown as ManifestJson, options).open();
}

export class PropertyDialog {
  readonly editor: ManifestEditor;
  readonly dialog: Dialog;

  private readonly _banner = div([], 'u2-manifest-banner');
  private readonly _status = span('', 'u2-manifest-status');
  private _saved = false;

  constructor(private readonly _schema: string, private readonly _handle: DG.DomainSchemaClient,
    manifest: ManifestJson, private readonly _options: EditPropertiesOptions) {
    this.dialog = Dialog.create(_options.title ?? PropertyDialog.titleOf(manifest, _options.table), {name: 'properties'});
    this.editor = this.dialog.runInScope(() => new ManifestEditor(null,
      {context: {mode: 'extend', storage: 'domain'}, baseline: manifest, treeShare: 0.45}));
    this.editor.root.classList.add('u2-manifest-editor-dialog');
    this._banner.hidden = true;
    const changes = computed(() => {
      this.editor.model.revision.value;
      return this.editor.model.changes();
    });
    this.dialog.add(divV([this._banner, this.editor.root, this._status], 'u2-manifest-properties'));
    this.dialog.effect(() => this._status.textContent = PropertyDialog.summary(changes.value));
    this.dialog.okDisabled(computed(() => changes.value.length === 0 ? 'No changes to save' : null));
    this.dialog.onOK(() => this._save(), {text: 'SAVE', primary: true});
  }

  private static readonly _extensible = new Map<string, Promise<Set<string>>>();

  /** The schema's tables the current user may add properties to — none without Extend; read once
   * per schema and session, as the opt-in changes only with a plugin deploy. */
  static extensible(schema: string): Promise<Set<string>> {
    let tables = PropertyDialog._extensible.get(schema);
    if (tables === undefined) {
      tables = PropertyDialog.canExtend(schema).then(async (can) => !can ? new Set<string>() :
        new Set(Object.entries((await grok.dapi.domains.schema(schema).manifest()).tables)
          .filter(([, t]) => t.extensible === true).map(([name]) => name)));
      tables.catch(() => PropertyDialog._extensible.delete(schema));
      PropertyDialog._extensible.set(schema, tables);
    }
    return tables;
  }

  /** Whether the current user holds Extend on the schema — what shows a "properties" entry point. */
  static async canExtend(schema: string): Promise<boolean> {
    const entity = (await grok.dapi.domains.schemas.list()).find((s) => s.name === schema);
    return entity !== undefined && entity.managedBy === 'package' &&
      await grok.dapi.permissions.check(entity, 'Extend');
  }

  static titleOf(manifest: ManifestJson, table?: string): string {
    const decl = table === undefined ? undefined : manifest.tables[table];
    const name = decl === undefined ? manifest.name : (decl['singularName'] as string | undefined) ?? decl.friendlyName ?? table!;
    return `${name.charAt(0).toUpperCase()}${name.slice(1)} properties`;
  }

  /** The footer line: what the save will do, by kind. */
  static summary(changes: ManifestChange[]): string {
    if (changes.length === 0)
      return 'No changes';
    const added = changes.filter((c) => / added( |$)/.test(c.text)).length;
    const removed = changes.filter((c) => c.removes).length;
    const changed = changes.length - added - removed;
    return [added > 0 ? plural(added, 'property added', 'properties added') : '',
      changed > 0 ? plural(changed, 'change', 'changes') : '',
      removed > 0 ? plural(removed, 'property removed', 'properties removed') : ''].filter((s) => s !== '').join(' · ');
  }

  open(): Promise<boolean> {
    return new Promise((resolve) => {
      this.dialog.onClose(() => resolve(this._saved));
      this.dialog.show({width: 980, height: 640});
      if (this._options.table !== undefined)
        void this.editor.tree.reveal(this._options.table);
    });
  }

  /** The dry run, the confirmation where values go, the apply; false keeps the dialog open. */
  private async _save(): Promise<boolean> {
    const {payload, changes} = this.editor.extendPlan();
    if (Object.keys(payload.extend).length === 0)
      return true;
    let plan: DG.DomainApplyPlan;
    try {
      plan = await this._handle.apply(payload, {dryRun: true});
    } catch (e) {
      const demanded = (e as {body?: {plan?: DG.DomainApplyPlan}} | null)?.body?.plan;
      if (DomainErrors.codeOf(e) !== 'destructive-confirmation-required' || demanded === undefined)
        return this._refused(e);
      plan = demanded;
    }
    let body: ExtendPayload = payload;
    if (plan.destructive) {
      if (!await PropertyDialog.confirmDrops(this.dialog, changes, plan))
        return false;
      body = {...payload, confirmDestructive: true};
    }
    try {
      await this._handle.apply(body);
    } catch (e) {
      return this._refused(e);
    }
    this._saved = true;
    const tables = Object.keys(payload.extend);
    grok.dapi.domains.invalidateUiCaches();
    notify.info(`Saved: ${PropertyDialog.summary(changes)}`);
    await this._options.onSaved?.(tables);
    return true;
  }

  /** A conflict reloads and replays; a refusal the server answered lands on the rows it names. */
  private async _refused(e: unknown): Promise<boolean> {
    if (DomainErrors.codeOf(e) === 'version-conflict') {
      const report = this.editor.rebase(await this._handle.manifest() as unknown as ManifestJson, null, null,
        {friendlyName: '', description: ''});
      this._say(PropertyDialog.reloaded(report));
      return false;
    }
    this.editor.diagnostics.value = BindingWizard.issues(e);
    throw new Error(DomainErrors.message(e));
  }

  private _say(text: string): void {
    this._banner.replaceChildren(span(text));
    this._banner.hidden = false;
  }

  static reloaded(report: RebaseReport): string {
    const lost = [...report.conflicts, ...report.dropped].map((op) => op.note ?? `${op.text} — no longer applies`);
    return 'Someone else changed these properties after you opened them, so they were reloaded with your edits ' +
      `on top. ${lost.length === 0 ? '' : `${lost.join('; ')}. `}Review and save again.`;
  }

  /** What a destructive save takes with it, one line per change; resolves to whether to go on. */
  static confirmDrops(owner: Dialog, changes: ManifestChange[], plan: DG.DomainApplyPlan): Promise<boolean> {
    const values = (c: ManifestChange): number => plan.drops.columns
      .find((d) => d.table === c.table && d.column === c.column)?.nonNullValues ?? 0;
    const lost = changes.filter((c) => c.removes).reduce((n, c) => n + values(c), 0);
    const line = (c: ManifestChange): HTMLElement => {
      const kind = c.removes ? 'remove' : / added( |$)/.test(c.text) ? 'add' : 'change';
      const n = c.removes ? values(c) : 0;
      return div([span(kind === 'add' ? '+' : kind === 'remove' ? '−' : '~', `u2-manifest-plan-sign u2-manifest-plan-${kind}`),
        span(c.text), n > 0 ? badge(`${plural(n, 'value', 'values')} deleted`, {variant: 'error'}) : span('')],
      'u2-manifest-plan-line');
    };
    return new Promise((resolve) => {
      const confirm = signal(lost === 0);
      const dialog = owner.runInScope(() => Dialog.create('Save these changes?', {name: 'properties-confirm'}));
      const box = dialog.runInScope(() => new BoolInput({label: lost === 0 ? 'Remove these properties' :
        `Delete ${plural(lost, 'stored value', 'stored values')} with them`, name: 'confirmDestructive',
      value: lost === 0, onChanged: (v) => confirm.value = v}));
      dialog.add(divV([span('Removing a property deletes its values from every row. This cannot be undone.'),
        divV(changes.map(line), 'u2-manifest-plan'), box.root], 'u2-manifest-confirm'));
      dialog.okDisabled(computed(() => confirm.value ? null : 'Confirm the deletion first'));
      let answer = false;
      dialog.onOK(() => answer = true, {text: 'SAVE', primary: true});
      dialog.onClose(() => resolve(answer));
      dialog.show({width: 520});
    });
  }
}

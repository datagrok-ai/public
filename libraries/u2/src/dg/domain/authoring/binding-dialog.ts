/* `domains.authoring.createBinding` — the "Create domain schema" dialog over an external database:
   **Connection** (the connection and its remote schema; the author's rights on it, the DML trio
   included; the draft is read as soon as both are picked, and what the server refuses is said
   right there), **Design** (the `ManifestEditor` over the draft) and **Review** (the manifest as
   it will be sent beside the findings; VALIDATE is the dry run, bound to the exact payload it
   checked; CREATE the real one, then the column restrictions, then the table grants). With access
   rows the dialog ends on **Created** — what was applied, what failed and why, RETRY ACCESS over
   the failed steps alone, OPEN — since a created schema is never reported as anything else. The
   server owns every rule — the dialog names its refusals. */
import * as grok from 'datagrok-api/grok';
import type * as DG from 'datagrok-api/dg';
import {computed, signal} from '../../../core/signals.js';
import {Control} from '../../../core/component.js';
import {div, divV, link, span} from '../../../core/elements.js';
import {plural} from '../../../core/text.js';
import {Wizard} from '../../../components/containers/wizard.js';
import {ChoiceInput} from '../../../components/inputs/choice-input.js';
import {Form} from '../../../components/forms/form.js';
import {badge} from '../../../components/display/badge.js';
import {notify} from '../../../components/display/notify.js';
import {DomainErrors} from '../errors.js';
import {route} from '../routes.js';
import {groupInput} from '../../inputs/group-input.js';
import {ManifestEditor} from './manifest-editor.js';
import type {ManifestPlan} from './manifest-editor.js';
import type {DraftEnvelope, ManifestDiagnostic, ManifestJson} from './manifest-model.js';

export interface CreateBindingOptions {
  connection?: DG.DataConnection;
  /** With `connection`, the dialog opens on Design; the Connection step stays a BACK away. */
  schema?: string;
  catalog?: string;
  /** Only this table starts included; the rest of the schema is drafted and offered excluded. */
  table?: string;
  /** The group names the access pickers find; every group and user the look-up answers otherwise. */
  groups?: string[];
}

/** How the access rows went after the create: each step by the line the Created state shows. */
export interface AccessOutcome {
  applied: string[];
  failed: {op: string, message: string}[];
}

export interface BindingResult {
  name: string;
  access: AccessOutcome;
}

/** The event `grok.events` carries once a schema is created and its access applied — for a JS
 * listener; the Dart Domains tree refreshes through its own entry point's return. */
export const SCHEMA_CREATED_EVENT = 'domain-schema-created';

/** Opens the dialog; resolves once it closes — to the created schema with its access outcome,
 * null where it was cancelled before the create. Refuses while domain databases are off. */
export function createBinding(options: CreateBindingOptions = {}): Promise<BindingResult | null> {
  if (grok.shell.settings.enableDomainDatabases !== true)
    return Promise.reject(new Error('Domain databases are a Beta feature — enable them in Settings > Beta'));
  return new BindingDialog(options).open();
}

/** The domain handle is a database to the platform, but not one to bind. */
const DOMAIN_SOURCE = 'Domain';

const DML = ['AddRows', 'ChangeValues', 'RemoveRows'];

/** One access step after the create: a column restriction, or one permission on one table. */
interface AccessOp {
  kind: 'restrict' | 'grant';
  table: string;
  label: string;
  run: () => Promise<unknown>;
}

export class BindingDialog extends Control {
  readonly wizard: Wizard;
  readonly connection: ChoiceInput;
  readonly schema: ChoiceInput;

  private readonly _options: CreateBindingOptions;
  private readonly _editor = signal<ManifestEditor | undefined>(undefined);
  private _connections: DG.DataConnection[] = [];
  private _takenNames: string[] = [];
  private readonly _draft = signal<DraftEnvelope | null>(null);
  /** What stands between the connection step and Design: the read in progress, or its refusal. */
  private readonly _reading = signal<string | null>(null);
  /** Why the editor could not be built over the draft. */
  private readonly _designProblem = signal<string | null>(null);
  /** The exact create payload the dry run passed; anything else is not validated. */
  private readonly _validated = signal<string | null>(null);
  private readonly _failed = signal(0);
  private readonly _facts = span('Only database connections are offered', 'u2-binding-facts');
  private readonly _access = span('', 'u2-binding-access');
  private readonly _designHost = div([], 'u2-binding-design');
  private readonly _json = document.createElement('pre');
  private readonly _issues = divV([], 'u2-binding-issues');
  private readonly _createdHost = divV([], 'u2-binding-created');
  private _built: DraftEnvelope | null = null;
  private _writeHint: string | undefined;
  private _describeGen = 0;
  private _readGen = 0;
  private _loaded: Promise<void> = Promise.resolve();
  private _described: Promise<void> = Promise.resolve();
  private _result: BindingResult | undefined;
  private _firstTable = '';
  private _pending: AccessOp[] = [];
  private _resolve: ((result: BindingResult | null) => void) | undefined;

  constructor(options: CreateBindingOptions = {}) {
    super();
    this._options = options;
    this.root.classList.add('u2-binding-dialog-host');
    this.root.dataset.u2 = 'binding-dialog';
    this.connection = this.runInScope(() => new ChoiceInput({label: 'Connection', name: 'connection', items: [],
      nullable: false, emptyText: 'Loading connections…'}));
    this.schema = this.runInScope(() => new ChoiceInput({label: 'Schema', name: 'schema', items: [],
      nullable: false, emptyText: 'Pick a connection first'}));
    this._json.className = 'u2-binding-json';
    // the connection step is built up front: it starts the reads, whichever step opens first
    const connection = this.runInScope(() => this._connectionStep());
    this.wizard = this.runInScope(() => new Wizard({
      start: options.connection !== undefined && options.schema !== undefined ? 'design' : undefined,
      steps: [
        {id: 'connection', title: 'Connection', content: connection,
          canProceed: () => this._draft.value !== null ? null :
            this._reading.value ?? 'Pick a connection and a schema'},
        {id: 'design', title: 'Design', content: () => this._designHost, canProceed: () => this._designGate()},
        {id: 'review', title: 'Review', content: () => this._reviewStep(), nextText: 'CREATE',
          onActivate: () => this._review(), actions: [{text: 'VALIDATE', run: () => this._validate()}],
          canProceed: () => this._reviewGate(), commit: () => this._create()},
        {id: 'created', title: 'Created', content: () => this._createdHost, done: true, actions: [
          {text: 'RETRY ACCESS', run: () => this._retry(), enabled: computed(() => this._failed.value > 0)},
          {text: 'OPEN', run: () => this._openAndClose()},
        ]},
      ],
      onFinish: () => this._settle(this._result ?? null),
      onCancel: () => this._settle(null),
    }));
    this.wizard.root.classList.add('u2-binding-dialog');
    this.root.append(this.wizard.root);
    // the editor is built once Design is open and the draft is in, whichever comes second
    this.effect(() => {
      this._draft.value;
      if (this.wizard.currentStep.value === 'design')
        void this._design();
    });
  }

  open(): Promise<BindingResult | null> {
    this.wizard.openInDialog('Create domain schema', {width: 960, height: 640});
    return new Promise((resolve) => this._resolve = resolve);
  }

  /** The editor, from the Design step on; rebuilt when the draft behind it changed. */
  get editor(): ManifestEditor | undefined {
    return this._editor.value;
  }

  /** The dry run's findings from a refusal: the manifest-validation errors by path, else the
   * refusal itself as one schema-wide finding. */
  static issues(e: unknown): ManifestDiagnostic[] {
    const errors = (e as {errors?: unknown} | null)?.errors;
    if (Array.isArray(errors) && errors.length > 0)
      return errors as ManifestDiagnostic[];
    return [{code: DomainErrors.codeOf(e) || 'error', message: DomainErrors.message(e)}];
  }

  private _settle(result: BindingResult | null): void {
    const resolve = this._resolve;
    this._resolve = undefined;
    this.wizard.dispose();
    this.dispose();
    resolve?.(result);
  }

  private _say(content: string | HTMLElement, error = false): void {
    this.wizard.status.replaceChildren(content);
    this.wizard.status.classList.toggle('u2-wizard-status-error', error);
  }

  private _picked(): DG.DataConnection | undefined {
    const id = this.connection.value.value;
    return this._connections.find((c) => c.id === id);
  }

  private _connectionStep(): HTMLElement {
    const form = new Form().addAll([this.connection, this.schema]);
    form.addElement(this._facts);
    form.addElement(this._access);
    const preset = this._options.connection;
    if (preset !== undefined)
      this._offerConnections([preset]);
    this._loaded = this._load(preset);
    // a new connection: its schemas, the author's rights on it; a new schema: the draft. Each
    // read is stamped, so the answer of a superseded one is dropped
    form.effect(() => {
      const conn = this._picked();
      const gen = ++this._describeGen;
      this.schema.setItems([]);
      this._access.textContent = '';
      if (conn !== undefined)
        this._described = this._describe(conn, gen);
    });
    form.effect(() => {
      const conn = this._picked();
      const schema = this.schema.value.value;
      const gen = ++this._readGen;
      this._draft.value = null;
      this._reading.value = null;
      if (conn !== undefined && schema !== null)
        void this._read(conn, schema, gen);
    });
    return divV([form.root], 'u2-binding-connection');
  }

  private _offerConnections(list: DG.DataConnection[]): void {
    const preset = this._options.connection;
    this._connections = preset !== undefined && !list.some((c) => c.id === preset.id) ? [preset, ...list] : list;
    this.connection.setItems(this._connections.map((c) =>
      ({value: c.id, label: `${c.friendlyName} (${c.dataSource})`})));
    if (preset !== undefined && this.connection.value.peek() === null)
      this.connection.value.value = preset.id;
  }

  private async _load(preset: DG.DataConnection | undefined): Promise<void> {
    let stage = 'Connections';
    try {
      const list = (await grok.dapi.connections.list({pageSize: 500}))
        .filter((c) => c.isDatabase && c.dataSource !== DOMAIN_SOURCE);
      list.sort((a, b) => a.friendlyName.localeCompare(b.friendlyName));
      this._offerConnections(list);
      stage = 'Registered schemas';
      this._takenNames = (await grok.dapi.domains.schemas.list({pageSize: 500})).map((s) => s.name);
    } catch (e) {
      if (stage === 'Connections' && preset === undefined)
        this.connection.setItems([]);
      this._facts.textContent = `${stage} could not be listed: ${DomainErrors.message(e)}`;
      this._facts.classList.add('u2-binding-problem');
    }
  }

  private async _describe(conn: DG.DataConnection, gen: number): Promise<void> {
    this._facts.className = 'u2-binding-facts';
    this._facts.textContent = `${conn.dataSource} · reading the schemas…`;
    this._writeHint = undefined;
    let stage = 'The schemas';
    try {
      const schemas = await grok.dapi.connections.getSchemas(conn, this._options.catalog ?? null);
      if (gen !== this._describeGen)
        return;
      this._facts.textContent = `${conn.dataSource} · ${plural(schemas.length, 'schema', 'schemas')}`;
      const preset = this._options.schema;
      const items = preset !== undefined && !schemas.includes(preset) ? [preset, ...schemas] : schemas;
      this.schema.setItems(items.map((s) => ({value: s, label: s})));
      if (preset !== undefined)
        this.schema.value.value = preset;
      stage = 'Your rights on the connection';
      const rights = ['GetSchema', 'Query', ...DML];
      const granted = await Promise.all(rights.map((r) => grok.dapi.permissions.check(conn, `DataConnection.${r}`)));
      if (gen !== this._describeGen)
        return;
      const missing = rights.filter((_, i) => !granted[i]);
      const read = missing.filter((r) => !DML.includes(r));
      const write = missing.filter((r) => DML.includes(r));
      this._writeHint = write.length === 0 ? undefined :
        `writes need ${write.join(', ')} on the connection, which you lack`;
      this._access.textContent = read.length > 0 ?
        `You lack ${read.join(' and ')} on this connection — the draft will be refused` :
        write.length === 0 ? 'You may introspect, query and write to this connection' :
          `You may introspect and query this connection; ${this._writeHint}`;
      this._access.classList.toggle('u2-binding-problem', read.length > 0);
    } catch (e) {
      if (gen !== this._describeGen)
        return;
      this._writeHint = 'your rights on the connection could not be read';
      this._facts.textContent = `${conn.dataSource} · ${stage} could not be read: ${DomainErrors.message(e)}`;
      this._facts.classList.add('u2-binding-problem');
    }
  }

  private async _read(conn: DG.DataConnection, schema: string, gen: number): Promise<void> {
    this._reading.value = 'Reading the schema…';
    try {
      const answer = await grok.dapi.domains.draft({connection: conn.nqName, schema, catalog: this._options.catalog});
      if (gen !== this._readGen)
        return;
      const draft: DraftEnvelope = {...answer, manifest: answer.manifest as ManifestJson};
      const tables = answer.inventory.tables;
      const bindable = tables.filter((t) => t.bindable).length;
      this._facts.textContent = `${conn.dataSource} · ${schema}: ${bindable} of ` +
        `${plural(tables.length, 'table', 'tables')} bindable`;
      this._draft.value = draft;
      this._reading.value = null;
    } catch (e) {
      if (gen !== this._readGen)
        return;
      this._reading.value = DomainErrors.message(e);
    }
  }

  /** The editor over the current draft, built on the first visit and again when the draft behind
   * it changed (BACK to the connection step, another schema). */
  private async _design(): Promise<void> {
    const draft = this._draft.peek();
    if (draft === null || draft === this._built)
      return;
    this._designProblem.value = null;
    try {
      await Promise.all([this._loaded, this._described]);
      if (this.scope.isDisposed || this._draft.peek() !== draft)
        return;
      this._editor.peek()?.dispose();
      this._built = draft;
      const wanted = this._options.groups;
      const editor = this.runInScope(() => new ManifestEditor(draft, {
        context: {mode: 'create', storage: 'external'}, takenNames: this._takenNames, writableDisabled: this._writeHint,
        principalPicker: (onPick) => groupInput({
          accept: wanted === undefined ? undefined : (g) => wanted.includes(g.friendlyName),
          onPick: (g, label) => onPick({id: g.id, label}),
        }).root,
      }));
      editor.model.setSchemaFriendlyName(`${this._picked()!.friendlyName} ${this.schema.value.peek()}`);
      const only = this._options.table;
      if (only !== undefined && editor.model.table(only) !== undefined) {
        editor.model.includeTables(false);
        editor.model.includeTable(only, true);
      }
      this._designHost.replaceChildren(editor.root);
      this._editor.value = editor;
    } catch (e) {
      this._designProblem.value = DomainErrors.message(e);
      this._say(this._designProblem.value, true);
    }
  }

  private _designGate(): string | null {
    const problem = this._designProblem.value;
    if (problem !== null)
      return problem;
    const editor = this._editor.value;
    if (editor === undefined)
      return this._reading.value ?? 'Reading the draft…';
    const name = editor.model.checkSchemaName(editor.model.name.value);
    if (name !== null)
      return `Name: ${name}`;
    return editor.model.tables.value.some((t) => t.included) ? null : 'Include at least one table';
  }

  private _reviewStep(): HTMLElement {
    return div([this._json, this._issues], 'u2-binding-review');
  }

  /** The create payload as one string — what a validation is bound to. */
  private static _payload(editor: ManifestEditor): string {
    const model = editor.model;
    return JSON.stringify({name: model.name.peek(), friendlyName: model.friendlyName.peek(), manifest: model.toJSON()});
  }

  private _isValidated(editor: ManifestEditor): boolean {
    const model = editor.model;
    model.revision.value;
    model.name.value;
    model.friendlyName.value;
    return this._validated.value === BindingDialog._payload(editor);
  }

  private _reviewGate(): string | null {
    const editor = this._editor.value;
    return editor !== undefined && this._isValidated(editor) ? null : 'Validate before creating';
  }

  private _review(): void {
    const editor = this.editor!;
    this._json.textContent = JSON.stringify(editor.model.toJSON(), null, 2);
    this._say(this._isValidated(editor) ? badge('Validated', {variant: 'success'}) : '');
    this._renderIssues();
  }

  private _renderIssues(): void {
    const editor = this.editor!;
    const issues = editor.diagnostics.peek();
    this._issues.replaceChildren(span(issues.length === 0 ? 'No findings' :
      plural(issues.length, 'finding', 'findings'), 'u2-binding-issues-title'));
    for (const issue of issues) {
      const row = div([link(issue.path ?? 'schema', () => {
        void editor.select(editor.model.resolvePath(issue.path));
        this.wizard.goTo('design');
      }), span(`: ${issue.message}`)], 'u2-binding-issue');
      this._issues.append(row);
    }
  }

  /** The dry run over the payload as it stands; an answer to a payload since changed is dropped. */
  private async _validate(): Promise<void> {
    const editor = this.editor!;
    const plan = editor.plan();
    const payload = BindingDialog._payload(editor);
    const stale = (): boolean => this.editor !== editor || BindingDialog._payload(editor) !== payload;
    this._say('Validating…');
    try {
      await grok.dapi.domains.createSchema(plan.name, {friendlyName: plan.friendlyName || undefined,
        manifest: plan.manifest, dryRun: true});
      if (stale())
        return this._say('Changed while validating — validate again');
      editor.diagnostics.value = [];
      this._validated.value = payload;
      this._say(badge('Validated', {variant: 'success'}));
    } catch (e) {
      if (stale())
        return this._say('Changed while validating — validate again');
      editor.diagnostics.value = BindingDialog.issues(e);
      this._validated.value = null;
      this._say(DomainErrors.message(e), true);
    }
    this._renderIssues();
  }

  /** The create, then the access rows; the schema exists once the create answered, and the
   * dialog reports it as created whatever the access rows say: without any, today's short way
   * (a toast, the app, the promise); with some, the Created step. */
  private async _create(): Promise<boolean> {
    const editor = this.editor!;
    const plan = editor.plan();
    this._say('Creating…');
    try {
      await grok.dapi.domains.createSchema(plan.name, {friendlyName: plan.friendlyName || undefined,
        manifest: plan.manifest});
    } catch (e) {
      editor.diagnostics.value = BindingDialog.issues(e);
      this._validated.value = null;
      this._renderIssues();
      this._say(DomainErrors.message(e), true);
      return false;
    }
    const ops = BindingDialog._accessOps(plan);
    this._result = {name: plan.name, access: {applied: [], failed: []}};
    this._firstTable = Object.keys(plan.manifest.tables)[0];
    this._say('Applying access…');
    await this._applyAccess(ops);
    grok.dapi.domains.invalidateUiCaches();
    grok.events.fireCustomEvent(SCHEMA_CREATED_EVENT, {name: plan.name});
    if (ops.length > 0) {
      this._renderCreated();
      return true;
    }
    notify.info(`Domain schema ${plan.name} created`);
    await this._openAndClose();
    return false;
  }

  private async _openAndClose(): Promise<void> {
    const result = this._result!;
    const first = this._firstTable;
    this._settle(result);
    try {
      const view = await route(`/domains/${result.name}/${first}`);
      if (view !== null)
        grok.shell.addView(view);
    } catch (e) {
      notify.error(`Domain schema ${result.name} could not be opened: ${DomainErrors.message(e)}`);
    }
  }

  private async _retry(): Promise<void> {
    this._say('Applying access…');
    await this._applyAccess(this._pending);
    this._renderCreated();
  }

  /** Restrictions first — a column meant for some is never readable by all in between — then the
   * grants, none on a table whose restriction failed. What failed stays for a retry. */
  private async _applyAccess(ops: AccessOp[]): Promise<void> {
    const access = this._result!.access;
    const failed: {op: AccessOp, message: string}[] = [];
    const broken = new Set<string>();
    const run = async (op: AccessOp): Promise<void> => {
      try {
        await op.run();
        access.applied.push(op.label);
      } catch (e) {
        failed.push({op, message: DomainErrors.message(e)});
        if (op.kind === 'restrict')
          broken.add(op.table);
      }
    };
    for (const op of ops.filter((o) => o.kind === 'restrict'))
      await run(op);
    for (const op of ops.filter((o) => o.kind === 'grant')) {
      if (broken.has(op.table))
        failed.push({op, message: `withheld — a column restriction of ${op.table} failed`});
      else
        await run(op);
    }
    this._pending = failed.map((f) => f.op);
    access.failed = failed.map((f) => ({op: f.op.label, message: f.message}));
    this._failed.value = failed.length;
  }

  private static _accessOps(plan: ManifestPlan): AccessOp[] {
    const table = (t: string): DG.DomainTableClient => grok.dapi.domains.table(`${plan.name}.${t}`);
    const ops: AccessOp[] = [];
    for (const r of plan.restrictions) {
      const who = r.groups.length === 0 ? 'nobody else' : r.groups.map((g) => g.label).join(', ');
      ops.push({kind: 'restrict', table: r.table, label: `${r.table}.${r.column} visible to ${who}`,
        run: async () => {
          if (r.groups.length === 0)
            await table(r.table).restrictColumn(r.column);
          for (const g of r.groups)
            await table(r.table).shareColumn(r.column, g.id, 'View');
        }});
    }
    for (const g of plan.grants) {
      const permissions = [...(g.view ? ['View'] : []), ...(g.edit ? ['Edit'] : []), ...(g.delete ? ['Delete'] : [])];
      for (const permission of permissions) {
        ops.push({kind: 'grant', table: g.table, label: `${g.table}: ${permission} for ${g.group.label}`,
          run: () => table(g.table).grant(g.group.id, permission)});
      }
    }
    return ops;
  }

  private _renderCreated(): void {
    const {name, access} = this._result!;
    const list = (title: string, lines: string[], cls: string): HTMLElement =>
      divV([span(title, 'u2-binding-created-title'), ...lines.map((l) => span(l))], `u2-binding-created-list ${cls}`);
    const content: HTMLElement[] = [div([badge('Created', {variant: 'success'}),
      span(`Domain schema ${name} is registered`)], 'u2-binding-created-head')];
    if (access.applied.length > 0)
      content.push(list('Access applied', access.applied, 'u2-binding-created-applied'));
    if (access.failed.length > 0) {
      content.push(list(`${plural(access.failed.length, 'access step', 'access steps')} failed — retry, or redo ` +
        'them from the schema\'s page', access.failed.map((f) => `${f.op}: ${f.message}`),
      'u2-binding-created-failed'));
    }
    this._createdHost.replaceChildren(...content);
    this._say(access.failed.length === 0 ? badge('Access applied', {variant: 'success'}) :
      `${plural(access.failed.length, 'access step', 'access steps')} failed`, access.failed.length > 0);
  }
}

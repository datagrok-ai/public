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
import {Permission} from 'datagrok-api/dg';
import {computed, signal} from '../../../core/signals.js';
import {Control} from '../../../core/component.js';
import {div, divV, link, span} from '../../../core/elements.js';
import {plural} from '../../../core/text.js';
import {Wizard} from '../../../components/containers/wizard.js';
import type {WizardOptions} from '../../../components/containers/wizard.js';
import {Section} from '../../../components/containers/section.js';
import {ChoiceInput} from '../../../components/inputs/choice-input.js';
import {Form} from '../../../components/forms/form.js';
import {badge} from '../../../components/display/badge.js';
import {notify} from '../../../components/display/notify.js';
import {DomainErrors} from '../errors.js';
import {route} from '../routes.js';
import {groupInput} from '../../inputs/group-input.js';
import {ManifestEditor} from './manifest-editor.js';
import type {ManifestPlan} from './manifest-editor.js';
import {ManifestModel} from './manifest-model.js';
import type {AccessPrincipal, DraftEnvelope, ManifestDiagnostic, ManifestJson} from './manifest-model.js';

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
  /** A create that got no answer, with a schema registered under the name: no access was applied. */
  unknown?: boolean;
}

/** The platform's Dart client learns of a create from the create itself; a schema found registered
 * after a lost answer is announced to it here as altered (read at call time like the designer's
 * globals). */
const api = globalThis as {grok_Dapi_Domains_SchemaAltered?: (dart: unknown, name: string) => void};

/** Opens the dialog; resolves once it closes — to the created schema with its access outcome,
 * null where it was cancelled before the create. Refuses while domain databases are off, and
 * without the privilege the create needs — by name, before anything is read. */
export async function createBinding(options: CreateBindingOptions = {}): Promise<BindingResult | null> {
  if (grok.shell.settings.enableDomainDatabases !== true)
    throw new Error('Domain databases are a Beta feature — enable them in Settings > Beta');
  if (!await grok.dapi.permissions.checkGlobal(Permission.CREATE_DOMAIN_SCHEMA))
    throw new Error(`Requires the ${Permission.CREATE_DOMAIN_SCHEMA} privilege`);
  return new BindingDialog(options).open();
}

/** The domain handle is a database to the platform, but not one to bind. */
const DOMAIN_SOURCE = 'Domain';

const DML = ['AddRows', 'ChangeValues', 'RemoveRows'];

const PAGE = 500;

/** How long a warehouse read (the schemas, the draft) may take before the dialog says so. */
const READ_TIMEOUT = 30_000;

/** One access step after the create: a column restriction or share to one group, or one
 * permission on one table. */
interface AccessOp {
  kind: 'restrict' | 'grant';
  table: string;
  label: string;
  group?: AccessPrincipal;
  run: () => Promise<unknown>;
}

/** A failed access step; `final` where retrying cannot help (the group is gone). */
interface AccessFailure {
  op: AccessOp;
  message: string;
  final: boolean;
}

/** What the create and the edit dialogs share: the wizard in a modal dialog, the status line,
 * the findings panel with its rows taken to Design, a validation bound to the exact payload it
 * checked (an edit since is not validated, and takes the findings off the list), the wording of a
 * failure the server did not answer, the read timeout, and OPEN over a registered schema. */
export abstract class BindingWizard<TResult> extends Control {
  static readTimeout = READ_TIMEOUT;

  wizard!: Wizard;

  protected readonly _editor = signal<ManifestEditor | undefined>(undefined);
  /** The exact payload the dry run passed; anything else is not validated. */
  protected readonly _validated = signal<string | null>(null);
  protected readonly _current = computed(() => {
    const editor = this._editor.value;
    return editor === undefined ? null : this._payloadOf(editor);
  });
  protected readonly _issues = divV([], 'u2-binding-issues');
  /** The payload the editor's diagnostics were found in. */
  protected _diagnosed: string | null = null;
  protected _validateGen = 0;
  private _resolve: ((result: TResult | null) => void) | undefined;

  constructor() {
    super();
    this.root.classList.add('u2-binding-dialog-host');
    // findings belong to the payload they were found in: an edit since takes them off the list,
    // the rows and the panel
    this.effect(() => {
      const payload = this._current.value;
      const editor = this._editor.peek();
      if (editor !== undefined && payload !== this._diagnosed && editor.diagnostics.peek().length > 0)
        editor.diagnostics.value = [];
    });
  }

  /** The editor, from the Design step on. */
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

  /** The payload a validation is bound to, as one string with sorted keys, read so that the
   * signals it depends on are tracked. */
  protected abstract _payloadOf(editor: ManifestEditor): string;

  protected _mount(options: WizardOptions): void {
    this.wizard = this.runInScope(() => new Wizard(options));
    this.wizard.root.classList.add('u2-binding-dialog');
    this.root.append(this.wizard.root);
  }

  protected _open(title: string): Promise<TResult | null> {
    this.wizard.openInDialog(title, {width: 960, height: 640});
    return new Promise((resolve) => this._resolve = resolve);
  }

  protected _settle(result: TResult | null): void {
    const resolve = this._resolve;
    this._resolve = undefined;
    this.wizard.dispose();
    this.dispose();
    resolve?.(result);
  }

  protected _say(content: string | HTMLElement, error = false): void {
    this.wizard.status.replaceChildren(content);
    this.wizard.status.classList.toggle('u2-wizard-status-error', error);
  }

  protected _isValidated(): boolean {
    const validated = this._validated.value;
    return validated !== null && validated === this._current.value;
  }

  protected _renderIssues(): void {
    const editor = this.editor!;
    const issues = editor.diagnostics.peek();
    this._issues.replaceChildren(span(issues.length > 0 ? plural(issues.length, 'finding', 'findings') :
      this._isValidated() ? 'No findings' : 'Not validated',
    'u2-binding-issues-title'));
    for (const issue of issues) {
      const row = div([link(issue.path ?? 'schema', () => {
        void editor.select(editor.model.resolvePath(issue.path));
        this.wizard.goTo('design');
      }), span(`: ${issue.message}`)], 'u2-binding-issue');
      this._issues.append(row);
    }
  }

  /** The `DomainSchema` entity of a binding (caption, description, author). The smart filter's
   * `name` is the caption and its free text searches the caption too, so the binding is found by
   * its physical schema, `ext_<name>`, which is exact. */
  protected static _entity(name: string): Promise<DG.DomainSchema | null> {
    return grok.dapi.domains.schemas.filter(`pgSchema = "ext_${name}"`).list()
      .then((list) => list.find((s) => s.name === name) ?? null);
  }

  /** The app over the table, once the dialog settled; a route no app answers is said. */
  protected static async _openApp(name: string, table: string): Promise<void> {
    try {
      const view = await route(`/domains/${name}/${table}`);
      if (view !== null)
        grok.shell.addView(view);
      else
        notify.error(`Domain schema ${name} is registered but could not be opened`);
    } catch (e) {
      notify.error(`Domain schema ${name} could not be opened: ${DomainErrors.message(e)}`);
    }
  }

  /** JSON with sorted keys: the same value always reads the same. */
  protected static _canonical(value: unknown): string {
    return JSON.stringify(value, (_, v) => v === null || typeof v !== 'object' || Array.isArray(v) ? v :
      Object.fromEntries(Object.keys(v).sort().map((k) => [k, v[k]])));
  }

  protected static _notFound(e: unknown): boolean {
    return (e as {status?: unknown} | null)?.status === 404 || DomainErrors.codeOf(e) === 'not-found';
  }

  /** A refusal the server answered — not a lost connection, a timeout or a 5xx. */
  protected static _answered(e: unknown): boolean {
    const status = (e as {status?: unknown} | null)?.status;
    return DomainErrors.codeOf(e) !== '' && !(typeof status === 'number' && status >= 500);
  }

  /** A warehouse read that outlasts {@link readTimeout} fails by name — a connector that never
   * answers (a catalog it cannot open) must not hold the dialog. */
  protected static _within<T>(read: Promise<T>, what: string): Promise<T> {
    let timer: ReturnType<typeof setTimeout>;
    const late = new Promise<never>((_, reject) => timer = setTimeout(() =>
      reject(new Error(`${what} did not answer within ${BindingWizard.readTimeout / 1000} s`)),
    BindingWizard.readTimeout));
    return Promise.race([read, late]).finally(() => clearTimeout(timer));
  }

  /** A failure for the status line: a domain call nothing answered (status 0, no code — the
   * transport's own words) says so around them, the bare transport failure of the platform's
   * client says so alone; anything else is the refusal as worded. */
  protected static _failure(e: unknown): string {
    const message = DomainErrors.message(e);
    if (/^XMLHttpRequest error\.?$/.test(message))
      return 'The connection to the server failed';
    return (e as {status?: unknown} | null)?.status === 0 && DomainErrors.codeOf(e) === '' ?
      `Could not reach the server (${message})` : message;
  }
}


export class BindingDialog extends BindingWizard<BindingResult> {
  readonly connection: ChoiceInput;
  readonly schema: ChoiceInput;

  private readonly _options: CreateBindingOptions;
  private _connections: DG.DataConnection[] = [];
  private readonly _draft = signal<DraftEnvelope | null>(null);
  /** What stands between the connection step and Design: the read in progress, or its refusal. */
  private readonly _reading = signal<string | null>(null);
  /** Why the editor could not be built over the draft. */
  private readonly _designProblem = signal<string | null>(null);
  private readonly _failed = signal(0);
  /** The create's outcome could not be established: the Created step says what is registered. */
  private readonly _unknown = signal(false);
  private readonly _openable = signal(true);
  private readonly _facts = span('Only database connections are offered', 'u2-binding-facts');
  private readonly _access = span('', 'u2-binding-access');
  private readonly _designHost = div([], 'u2-binding-design');
  private readonly _json = document.createElement('pre');
  private readonly _planned: Section;
  private readonly _createdHost = divV([], 'u2-binding-created');
  private _built: DraftEnvelope | null = null;
  /** The connection and the remote schema the built draft was read over. */
  private _builtKey: string | null = null;
  /** The name of a create that got no answer: a later "name taken" may be that create's own. */
  private _unanswered: string | null = null;
  private _writeHint: string | undefined;
  /** A failed draft read is on the status line, to be cleared by the next read. */
  private _readProblem = false;
  private _describeGen = 0;
  private _readGen = 0;
  private _loaded: Promise<void> = Promise.resolve();
  private _described: Promise<void> = Promise.resolve();
  private _result: BindingResult | undefined;
  private _firstTable = '';
  private _pending: AccessOp[] = [];
  private readonly _final: AccessFailure[] = [];

  constructor(options: CreateBindingOptions = {}) {
    super();
    this._options = options;
    this.root.dataset.u2 = 'binding-dialog';
    this.connection = this.runInScope(() => new ChoiceInput({label: 'Connection', name: 'connection', items: [],
      nullable: false, emptyText: 'Loading connections…'}));
    this.schema = this.runInScope(() => new ChoiceInput({label: 'Schema', name: 'schema', items: [],
      nullable: false, emptyText: 'Pick a connection first'}));
    this._json.className = 'u2-binding-json';
    this._planned = this.runInScope(() => new Section({title: 'Planned access', collapsible: false}));
    this._planned.root.classList.add('u2-binding-planned');
    // the connection step is built up front: it starts the reads, whichever step opens first
    const connection = this.runInScope(() => this._connectionStep());
    this._mount({
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
          {text: 'RETRY ACCESS', run: () => this._retry(), enabled: computed(() => this._failed.value > 0),
            visible: computed(() => !this._unknown.value)},
          {text: 'OPEN', run: () => this._openAndClose(), visible: this._openable},
        ]},
      ],
      onFinish: () => this._settle(this._result ?? null),
      onCancel: () => this._settle(null),
    });
    // the editor is built once Design is open and the draft is in, whichever comes second
    this.effect(() => {
      this._draft.value;
      if (this.wizard.currentStep.value === 'design')
        void this._design();
    });
  }

  open(): Promise<BindingResult | null> {
    return this._open('Create domain schema');
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
      if (this._validated.peek() !== null) {
        this._validated.value = null;
        this._say('');
      }
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

  /** Every connection, page by page: the list is not capped. */
  private async _load(preset: DG.DataConnection | undefined): Promise<void> {
    try {
      const list: DG.DataConnection[] = [];
      for (let page = 1; ; page++) {
        const batch = await grok.dapi.connections.list({pageSize: PAGE, pageNumber: page, order: 'id'});
        list.push(...batch.filter((c) => c.isDatabase && c.dataSource !== DOMAIN_SOURCE));
        if (batch.length < PAGE)
          break;
      }
      list.sort((a, b) => a.friendlyName.localeCompare(b.friendlyName));
      this._offerConnections(list);
    } catch (e) {
      if (preset === undefined)
        this.connection.setItems([]);
      this._facts.textContent = `Connections could not be listed: ${BindingDialog._failure(e)}`;
      this._facts.classList.add('u2-binding-problem');
    }
  }

  /** The catalog preset belongs to the preset connection; another one picked is read whole. */
  private _catalogOf(conn: DG.DataConnection): string | undefined {
    return conn.id === this._options.connection?.id ? this._options.catalog : undefined;
  }

  private async _describe(conn: DG.DataConnection, gen: number): Promise<void> {
    this._facts.className = 'u2-binding-facts';
    this._facts.textContent = `${conn.dataSource} · reading the schemas…`;
    this._writeHint = undefined;
    let stage = 'The schemas';
    try {
      const schemas = [...await BindingDialog._within(grok.dapi.connections.getSchemas(conn, this._catalogOf(conn) ?? null),
        'The connection')].sort((a, b) => a.toLowerCase().localeCompare(b.toLowerCase()));
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
      this._facts.textContent = `${conn.dataSource} · ${stage} could not be read: ${BindingDialog._failure(e)}`;
      this._facts.classList.add('u2-binding-problem');
      // a preset schema is still worth a draft: the read names its own failure, the picker stays
      const preset = this._options.schema;
      if (stage === 'The schemas' && preset !== undefined && this.schema.value.peek() === null) {
        this.schema.setItems([{value: preset, label: preset}]);
        this.schema.value.value = preset;
      }
    }
  }

  private async _read(conn: DG.DataConnection, schema: string, gen: number): Promise<void> {
    this._reading.value = 'Reading the schema…';
    if (this._readProblem) {
      this._readProblem = false;
      this._say('');
    }
    try {
      const answer = await BindingDialog._within(
        grok.dapi.domains.draft({connection: conn.nqName, schema, catalog: this._catalogOf(conn)}), 'The draft');
      if (gen !== this._readGen)
        return;
      const draft: DraftEnvelope = {...answer, manifest: answer.manifest as ManifestJson};
      const tables = answer.inventory.tables;
      const bindable = tables.filter((t) => t.bindable).length;
      this._facts.textContent = `${conn.dataSource} · ${schema}: ${bindable} of ` +
        `${plural(tables.length, 'table', 'tables')} bindable`;
      // the pick the editor was built over, picked again: its edits stand
      this._draft.value = BindingDialog._key(conn, schema) === this._builtKey ? this._built : draft;
      this._reading.value = null;
    } catch (e) {
      if (gen !== this._readGen)
        return;
      this._reading.value = BindingDialog._failure(e);
      // Design opened on the preset shows the gate's reason under NEXT; the failure belongs on the line too
      if (this.wizard.currentStep.peek() === 'design') {
        this._readProblem = true;
        this._say(this._reading.value, true);
      }
    }
  }

  /** The editor over the current draft, built on the first visit and again when the draft behind
   * it changed (another connection or schema; the same one picked again keeps it). */
  private async _design(): Promise<void> {
    const draft = this._draft.peek();
    if (draft === null || draft === this._built)
      return;
    this._designProblem.value = null;
    try {
      await Promise.all([this._loaded, this._described]);
      if (this.scope.isDisposed || this._draft.peek() !== draft)
        return;
      const conn = this._picked()!;
      const schema = this.schema.value.peek()!;
      const wanted = this._options.groups;
      const editor = this.runInScope(() => new ManifestEditor(draft, {
        context: {mode: 'create', storage: 'external'}, writableDisabled: this._writeHint,
        principalPicker: (onPick) => groupInput({
          accept: wanted === undefined ? undefined : (g) => wanted.includes(g.friendlyName),
          onPick: (g, label) => onPick({id: g.id, label}),
        }).root,
      }));
      const label = `${conn.friendlyName} ${schema}`;
      const unchecked = await this._clearName(editor.model, label);
      if (this.scope.isDisposed || this._draft.peek() !== draft) {
        editor.dispose();
        return;
      }
      editor.model.setSchemaFriendlyName(editor.model.proposeFriendlyName(label));
      const only = this._options.table;
      // the remote name as the catalog reports it, whatever case the caller had it in
      const table = only === undefined ? undefined :
        editor.model.tables.peek().find((t) => t.remote.toLowerCase() === only.toLowerCase());
      if (table !== undefined) {
        editor.model.includeTables(false);
        editor.model.includeTable(table.remote, true);
      }
      this._diagnosed = this._payloadOf(editor);
      const reset = this._editor.peek();
      reset?.dispose();
      this._designHost.replaceChildren(editor.root);
      this._built = draft;
      this._builtKey = BindingDialog._key(conn, schema);
      this._editor.value = editor;
      const notes: string[] = [];
      if (reset !== undefined)
        notes.push(`The design was reset: ${conn.friendlyName} · ${schema} is a new draft`);
      if (only !== undefined && table === undefined)
        notes.push(`Table ${only} is not in ${schema} — every table starts included`);
      if (unchecked !== null)
        notes.push(unchecked);
      if (notes.length > 0)
        this._say(notes.join('; '), unchecked !== null);
    } catch (e) {
      this._designProblem.value = BindingDialog._failure(e);
      this._say(this._designProblem.value, true);
    }
  }

  /** The registry is not listed — it may hold thousands of schemas: the proposed identifier is
   * probed by name and stepped past every registered one. Answers why the probe stopped short,
   * null when it settled: an unanswered probe leaves the name as proposed (a taken one is
   * refused at CREATE all the same). */
  private async _clearName(model: ManifestModel, label: string): Promise<string | null> {
    for (;;) {
      const name = model.proposeName(label);
      if (name === '')
        return null;
      try {
        await grok.dapi.domains.schema(name).manifest();
      } catch (e) {
        return BindingDialog._notFound(e) ? null : `Registered schemas could not be checked: ${BindingDialog._failure(e)}`;
      }
      model.markTaken(name);
    }
  }

  private _designGate(): string | null {
    const problem = this._designProblem.value;
    if (problem !== null)
      return problem;
    const editor = this._editor.value;
    const draft = this._draft.value;
    if (editor === undefined || draft === null || draft !== this._built)
      return this._reading.value ?? 'Reading the draft…';
    const name = editor.model.checkSchemaName(editor.model.name.value);
    if (name !== null)
      return `Name: ${name}`;
    return editor.model.tables.value.some((t) => t.included) ? null : 'Include at least one table';
  }

  private _reviewStep(): HTMLElement {
    return div([this._json, divV([this._planned.root, this._issues], 'u2-binding-review-side')], 'u2-binding-review');
  }

  /** The create envelope: the name, the friendly name and the manifest — the access rows are
   * applied afterwards and are no part of it. */
  protected _payloadOf(editor: ManifestEditor): string {
    const model = editor.model;
    model.revision.value;
    return BindingDialog._canonical({name: model.name.value, friendlyName: model.friendlyName.value,
      manifest: model.toJSON()});
  }

  private static _key(conn: DG.DataConnection, schema: string): string {
    return `${conn.id}\u0000${schema}`;
  }

  /** The validation holds for the editor over the draft picked now, nothing older. */
  private _reviewGate(): string | null {
    const draft = this._draft.value;
    return draft !== null && draft === this._built && this._isValidated() ? null : 'Validate before creating';
  }

  /** Rebuilt on every activation — the step's content is built once, and Design may have changed
   * since. */
  private _review(): void {
    const editor = this.editor!;
    this._json.textContent = JSON.stringify(editor.model.toJSON(), null, 2);
    this._say(this._isValidated() ? badge('Validated', {variant: 'success'}) : '');
    this._renderPlanned();
    this._renderIssues();
  }

  /** What CREATE applies after the create, in the order it applies it, by the labels the Created
   * step reports — the dry run checks the manifest alone. */
  private _renderPlanned(): void {
    const ops = BindingDialog._accessOps(this.editor!.plan());
    this._planned.body.replaceChildren(
      span(ops.length === 0 ? 'No additional table or column grants planned' :
        'Applied after CREATE, in this order. Not covered by Validate.',
      'u2-binding-planned-note'),
      ...ops.map((op) => span(op.label, 'u2-binding-planned-op')));
  }

  /** The dry run over the payload as it stands; an answer to a payload since changed, or to an
   * older run, is dropped. A failure the server did not answer is no finding. */
  private async _validate(): Promise<void> {
    const editor = this.editor!;
    const plan = editor.plan();
    const payload = this._payloadOf(editor);
    const gen = ++this._validateGen;
    this._say('Validating…');
    let refusal: unknown = null;
    try {
      await grok.dapi.domains.createSchema(plan.name, {friendlyName: plan.friendlyName || undefined,
        manifest: plan.manifest, dryRun: true});
    } catch (e) {
      refusal = e;
    }
    if (gen !== this._validateGen)
      return;
    if (this.editor !== editor || this._payloadOf(editor) !== payload)
      return this._say('Changed while validating — validate again');
    if (refusal !== null && !BindingDialog._answered(refusal))
      return this._say(BindingDialog._failure(refusal), true);
    this._diagnosed = payload;
    editor.diagnostics.value = refusal === null ? [] : BindingDialog.issues(refusal);
    this._validated.value = refusal === null ? payload : null;
    this._say(refusal === null ? badge('Validated', {variant: 'success'}) : DomainErrors.message(refusal),
      refusal !== null);
    this._renderIssues();
  }

  /** The create, then the access rows; the schema exists once the create answered, and the
   * dialog reports it as created whatever the access rows say: without any, today's short way
   * (a toast, the app, the promise); with some, the Created step. */
  private async _create(): Promise<boolean> {
    const plan = this.editor!.plan();
    this._say('Creating…');
    let missing: string[];
    try {
      missing = await BindingDialog._missingGroups(plan);
    } catch (e) {
      this._say(`The access groups could not be checked: ${BindingDialog._failure(e)}`, true);
      return false;
    }
    if (missing.length > 0) {
      const one = missing.length === 1;
      this._say(`${one ? `Group ${missing[0]} no longer exists` : `Groups ${missing.join(', ')} no longer exist`} — ` +
        `take ${one ? 'it' : 'them'} out of the access rows`, true);
      return false;
    }
    try {
      await grok.dapi.domains.createSchema(plan.name, {friendlyName: plan.friendlyName || undefined,
        manifest: plan.manifest});
    } catch (e) {
      return await this._lost(plan, e);
    }
    this._unanswered = null;
    const ops = BindingDialog._accessOps(plan);
    this._result = {name: plan.name, access: {applied: [], failed: []}};
    this._firstTable = Object.keys(plan.manifest.tables)[0];
    this._say('Applying access…');
    await this._applyAccess(ops);
    grok.dapi.domains.invalidateUiCaches();
    if (ops.length > 0) {
      this._renderCreated();
      return true;
    }
    notify.info(`Domain schema ${plan.name} created`);
    await this._openAndClose();
    return false;
  }

  /** Whether a create that failed ends the dialog with the outcome unknown: its answer lost, or a
   * create of ours after such a loss refused as taken, and a schema registered under the name —
   * whatever it holds, nothing proves it this create's own, so no access is applied. A refusal
   * lands on the findings; a failure the registry does not explain on the status line, CREATE
   * still offered. */
  private async _lost(plan: ManifestPlan, e: unknown): Promise<boolean> {
    const name = plan.name;
    const ours = DomainErrors.codeOf(e) === 'schema-name-taken' && this._unanswered === name;
    if (!BindingDialog._answered(e) || ours) {
      this._unanswered = name;
      this._say('No answer to the create — reading the registry…');
      let registered: DG.DomainRegisteredManifest | null = null;
      try {
        registered = await grok.dapi.domains.schema(name).manifest();
      } catch (x) {
        if (!BindingDialog._notFound(x)) {
          this._say(`${BindingDialog._failure(e)} — the registry could not be read (${BindingDialog._failure(x)}); CREATE again`, true);
          return false;
        }
      }
      if (registered !== null) {
        this._outcomeUnknown(plan, e, registered);
        return true;
      }
      if (!ours) {
        this._say(`${BindingDialog._failure(e)} — ${name} is not registered; CREATE again`, true);
        return false;
      }
      // taken, yet gone by the time the registry was read: an ordinary refusal
      this._unanswered = null;
    }
    const editor = this.editor!;
    this._diagnosed = this._payloadOf(editor);
    editor.diagnostics.value = BindingDialog.issues(e);
    this._validated.value = null;
    this._renderIssues();
    this._say(DomainErrors.message(e), true);
    return false;
  }

  /** No access is sent; the Created step says what is registered under the name. */
  private _outcomeUnknown(plan: ManifestPlan, e: unknown, registered: DG.DomainRegisteredManifest): void {
    const tables = Object.keys(registered.tables ?? {});
    const storage = registered.storage as {connection?: string, schema?: string} | undefined;
    this._unanswered = null;
    grok.dapi.domains.invalidateUiCaches();
    api.grok_Dapi_Domains_SchemaAltered?.(grok.dapi.domains.dart, plan.name);
    this._unknown.value = true;
    this._openable.value = tables.length > 0;
    this._result = {name: plan.name, access: {applied: [], failed: []}, unknown: true};
    this._firstTable = tables[0] ?? '';
    const facts = [
      `No answer to the create (${BindingDialog._failure(e)}) — whether it landed is unknown`,
      `${plan.name} is registered at version ${registered.version}` +
        `${storage?.connection === undefined ? '' : ` over ${storage.connection} · ${storage.schema ?? ''}`}, ` +
        `${plural(tables.length, 'table', 'tables')}${tables.length > 0 ? ` (${tables.join(', ')})` : ''}`,
      'Access was not applied — once you know the binding is yours, use Edit binding… to set it',
    ];
    this._createdHost.replaceChildren(
      div([badge('Outcome unknown', {variant: 'warning'}), span(`A schema ${plan.name} is registered`)], 'u2-binding-created-head'),
      divV(facts.map((f) => span(f)), 'u2-binding-created-list u2-binding-created-unknown'));
    this._say('Outcome unknown — no access was applied', true);
  }

  /** The groups of the access rows the platform no longer knows, by label, found in one look-up. */
  private static async _missingGroups(plan: ManifestPlan): Promise<string[]> {
    const labels = new Map<string, string>();
    for (const g of plan.grants)
      labels.set(g.group.id, g.group.label);
    for (const r of plan.restrictions) {
      for (const g of r.groups)
        labels.set(g.id, g.label);
    }
    if (labels.size === 0)
      return [];
    const ids = [...labels.keys()];
    const found = await grok.dapi.groups.list({filter: `id in (${ids.map((id) => `"${id}"`).join(', ')})`,
      pageSize: ids.length});
    const known = new Set(found.map((g) => g.id));
    return ids.filter((id) => !known.has(id)).map((id) => labels.get(id)!);
  }

  private async _openAndClose(): Promise<void> {
    const result = this._result!;
    const first = this._firstTable;
    this._settle(result);
    await BindingDialog._openApp(result.name, first);
  }

  private async _retry(): Promise<void> {
    this._say('Applying access…');
    await this._applyAccess(this._pending);
    this._renderCreated();
  }

  /** Restrictions first — a column meant for some is never readable by all in between — then the
   * grants, none on a table whose restriction failed. What failed stays for a retry, unless its
   * group is gone: no retry brings that back. */
  private async _applyAccess(ops: AccessOp[]): Promise<void> {
    const access = this._result!.access;
    const failed: AccessFailure[] = [];
    // a table whose column restriction can never be applied keeps its grants withheld for good
    const dead = new Set(this._final.filter((f) => f.op.kind === 'restrict').map((f) => f.op.table));
    const broken = new Set(dead);
    const run = async (op: AccessOp): Promise<void> => {
      try {
        await op.run();
        access.applied.push(op.label);
      } catch (e) {
        const message = DomainErrors.message(e);
        const gone = op.group !== undefined && /Unknown group/i.test(message);
        failed.push({op, message: gone ? `group ${op.group!.label} no longer exists` : message, final: gone});
        if (op.kind === 'restrict')
          broken.add(op.table);
        if (op.kind === 'restrict' && gone)
          dead.add(op.table);
      }
    };
    for (const op of ops.filter((o) => o.kind === 'restrict'))
      await run(op);
    for (const op of ops.filter((o) => o.kind === 'grant')) {
      if (broken.has(op.table)) {
        failed.push({op, message: `withheld — a column restriction of ${op.table} failed`,
          final: dead.has(op.table)});
      } else
        await run(op);
    }
    const open = failed.filter((f) => !f.final);
    this._final.push(...failed.filter((f) => f.final));
    this._pending = open.map((f) => f.op);
    access.failed = [...this._final, ...open].map((f) => ({op: f.op.label, message: f.message}));
    this._failed.value = open.length;
  }

  /** One step per column, group and permission, so a retry repeats no share that went through; a
   * visibility group gets View, and Edit too where the table is editable, or it could not write it.
   * The labels are what Review previews and Created reports. */
  private static _accessOps(plan: ManifestPlan): AccessOp[] {
    const table = (t: string): DG.DomainTableClient => grok.dapi.domains.table(`${plan.name}.${t}`);
    const ops: AccessOp[] = [];
    for (const r of plan.restrictions) {
      const column = `${r.table}.${r.column}`;
      if (r.groups.length === 0) {
        ops.push({kind: 'restrict', table: r.table, label: `${column} visible to nobody else`,
          run: () => table(r.table).restrictColumn(r.column)});
      }
      // Edit on the column writes nothing without Edit on the table, which may come through any group
      const edit = plan.grants.some((x) => x.table === r.table && x.edit);
      for (const g of r.groups) {
        for (const permission of edit ? ['View', 'Edit'] : ['View']) {
          ops.push({kind: 'restrict', table: r.table, group: g,
            label: permission === 'Edit' ? `${column} editable by ${g.label} (with table Edit)` :
              `${column} visible to ${g.label}`,
            run: () => table(r.table).shareColumn(r.column, g.id, permission)});
        }
      }
    }
    for (const g of plan.grants) {
      const permissions = [...(g.view ? ['View'] : []), ...(g.edit ? ['Edit'] : []), ...(g.delete ? ['Delete'] : [])];
      for (const permission of permissions) {
        ops.push({kind: 'grant', table: g.table, group: g.group,
          label: `${g.table}: ${permission} for ${g.group.label}`,
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
      content.push(list(`${plural(access.failed.length, 'access step', 'access steps')} failed — ` +
        `${this._failed.peek() > 0 ? 'retry, or redo' : 'redo'} them from the schema's page`,
      access.failed.map((f) => `${f.op}: ${f.message}`),
      'u2-binding-created-failed'));
    }
    this._createdHost.replaceChildren(...content);
    this._say(access.failed.length === 0 ? badge('Access applied', {variant: 'success'}) :
      `${plural(access.failed.length, 'access step', 'access steps')} failed`, access.failed.length > 0);
  }
}

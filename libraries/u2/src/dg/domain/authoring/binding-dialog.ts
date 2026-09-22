/* `domains.authoring.createBinding` — the "Create domain schema" dialog over an external database:
   **Connection** (the connection and its remote schema; the author's rights on it; the draft is
   read as soon as both are picked, and what the server refuses is said right there), **Design**
   (the `ManifestEditor` over the draft) and **Review** (the manifest as it will be sent beside
   the findings; VALIDATE is the dry run, CREATE the real one, then the access rows, then the app
   over the first table). The server owns every rule — the dialog names its refusals. */
import * as grok from 'datagrok-api/grok';
import type * as DG from 'datagrok-api/dg';
import {signal} from '../../../core/signals.js';
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
import {ManifestEditor} from './manifest-editor.js';
import type {ManifestPlan} from './manifest-editor.js';
import type {DraftEnvelope, ManifestDiagnostic, ManifestJson} from './manifest-model.js';

export interface CreateBindingOptions {
  connection?: DG.DataConnection;
  schema?: string;
  catalog?: string;
  /** Only this table starts included; the rest of the schema is drafted and offered excluded. */
  table?: string;
  /** The group names the access pickers offer; every non-personal group otherwise. */
  groups?: string[];
}

/** Opens the dialog; resolves to the created schema's name, or null where it was cancelled. */
export function createBinding(options: CreateBindingOptions = {}): Promise<string | null> {
  return new BindingDialog(options).open();
}

/** The platform's file shares, secret stores and the domain handle: what `DataSourceType.isDatabase`
 * (connection_info.dart) rules out, kept here since the js-api has no such helper. */
const NON_DATABASE_SOURCES = new Set(['AzureBlob', 'Dropbox', 'Files', 'GitHub', 'GoogleCloud', 'S3', 'CoreWeave',
  'Git', 'SharePoint', 'EFS', 'AWS', 'GCP', 'Domain']);

export class BindingDialog extends Control {
  readonly wizard: Wizard;
  readonly connection: ChoiceInput;
  readonly schema: ChoiceInput;

  private readonly _options: CreateBindingOptions;
  private readonly _editor = signal<ManifestEditor | undefined>(undefined);
  private _connections: DG.DataConnection[] = [];
  /** Group name → id; null where two groups carry the name. */
  private readonly _groups = new Map<string, string | null>();
  private readonly _draft = signal<DraftEnvelope | null>(null);
  /** What stands between the connection step and Design: the read in progress, or its refusal. */
  private readonly _reading = signal<string | null>(null);
  /** Why the editor could not be built over the draft. */
  private readonly _designProblem = signal<string | null>(null);
  private readonly _validated = signal(false);
  private readonly _facts = span('Only database connections are offered', 'u2-binding-facts');
  private readonly _access = span('', 'u2-binding-access');
  private readonly _designHost = div([], 'u2-binding-design');
  private readonly _json = document.createElement('pre');
  private readonly _issues = divV([], 'u2-binding-issues');
  private _built: DraftEnvelope | null = null;
  private _describeGen = 0;
  private _readGen = 0;
  private _loaded: Promise<void> = Promise.resolve();
  private _resolve: ((name: string | null) => void) | undefined;

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
    this.wizard = this.runInScope(() => new Wizard({
      steps: [
        {id: 'connection', title: 'Connection', content: () => this._connectionStep(),
          canProceed: () => this._draft.value !== null ? null :
            this._reading.value ?? 'Pick a connection and a schema'},
        {id: 'design', title: 'Design', content: () => this._designHost,
          onActivate: () => void this._design(), canProceed: () => this._designGate()},
        {id: 'review', title: 'Review', content: () => this._reviewStep(), finishText: 'CREATE',
          onActivate: () => this._review(), actions: [{text: 'VALIDATE', run: () => this._validate()}],
          canProceed: () => this._validated.value ? null : 'Validate before creating'},
      ],
      onFinish: () => this._create(),
      onCancel: () => this._settle(null),
    }));
    this.wizard.root.classList.add('u2-binding-dialog');
    this.root.append(this.wizard.root);
  }

  open(): Promise<string | null> {
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

  private _settle(name: string | null): void {
    const resolve = this._resolve;
    this._resolve = undefined;
    this.wizard.dispose();
    this.dispose();
    resolve?.(name);
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
        void this._describe(conn, gen);
    });
    form.effect(() => {
      const conn = this._picked();
      const schema = this.schema.value.value;
      const gen = ++this._readGen;
      this._draft.value = null;
      this._validated.value = false;
      this._reading.value = null;
      if (conn !== undefined && schema !== null)
        void this._read(conn, schema, gen);
    });
    return divV([form.root], 'u2-binding-connection');
  }

  private _offerConnections(list: DG.DataConnection[]): void {
    const preset = this._options.connection;
    this._connections = preset !== undefined && !list.some((c) => c.id === preset.id) ? [preset, ...list] : list;
    this.connection.setItems(this._connections.map((c) => ({value: c.id, label: `${c.friendlyName} (${c.dataSource})`})));
    if (preset !== undefined && this.connection.value.peek() === null)
      this.connection.value.value = preset.id;
  }

  private async _load(preset: DG.DataConnection | undefined): Promise<void> {
    let stage = 'Connections';
    try {
      const list = (await grok.dapi.connections.list({pageSize: 500}))
        .filter((c) => !NON_DATABASE_SOURCES.has(c.dataSource));
      list.sort((a, b) => a.friendlyName.localeCompare(b.friendlyName));
      this._offerConnections(list);
      stage = 'Groups';
      for (const g of await grok.dapi.groups.filter('personal = false').list({pageSize: 500}))
        this._groups.set(g.friendlyName, this._groups.has(g.friendlyName) ? null : g.id);
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
      const [getSchema, query] = await Promise.all([grok.dapi.permissions.check(conn, 'DataConnection.GetSchema'),
        grok.dapi.permissions.check(conn, 'DataConnection.Query')]);
      if (gen !== this._describeGen)
        return;
      const missing = [...(getSchema ? [] : ['GetSchema']), ...(query ? [] : ['Query'])];
      this._access.textContent = missing.length === 0 ? 'You may introspect and query this connection' :
        `You lack ${missing.join(' and ')} on this connection — the draft will be refused`;
      this._access.classList.toggle('u2-binding-problem', missing.length > 0);
    } catch (e) {
      if (gen !== this._describeGen)
        return;
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
      await this._loaded;
      if (this.scope.isDisposed || this._draft.peek() !== draft)
        return;
      this._editor.peek()?.dispose();
      this._built = draft;
      this._validated.value = false;
      const offered = [...this._groups.keys()];
      const editor = this.runInScope(() => new ManifestEditor(draft, {
        context: {mode: 'create', storage: 'external'},
        groups: this._options.groups === undefined ? offered : offered.filter((g) => this._options.groups!.includes(g)),
      }));
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
      return 'Reading the draft…';
    const name = editor.model.checkSchemaName(editor.model.name.value);
    if (name !== null)
      return `Name: ${name}`;
    return editor.model.tables.value.some((t) => t.included) ? null : 'Include at least one table';
  }

  private _reviewStep(): HTMLElement {
    return div([this._json, this._issues], 'u2-binding-review');
  }

  private _review(): void {
    const editor = this.editor!;
    this._json.textContent = JSON.stringify(editor.model.toJSON(), null, 2);
    this._renderIssues();
  }

  private _renderIssues(): void {
    const editor = this.editor!;
    const issues = editor.diagnostics.peek();
    this._issues.replaceChildren(span(issues.length === 0 ? 'No findings' : plural(issues.length, 'finding', 'findings'),
      'u2-binding-issues-title'));
    for (const issue of issues) {
      const row = div([link(issue.path ?? 'schema', () => {
        void editor.select(editor.model.resolvePath(issue.path));
        this.wizard.goTo('design');
      }), span(`: ${issue.message}`)], 'u2-binding-issue');
      this._issues.append(row);
    }
  }

  private async _validate(): Promise<void> {
    const editor = this.editor!;
    const plan = editor.plan();
    this._say('Validating…');
    try {
      await grok.dapi.domains.createSchema(plan.name, {friendlyName: plan.friendlyName || undefined,
        manifest: plan.manifest, dryRun: true});
      editor.diagnostics.value = [];
      this._validated.value = true;
      this._say(badge('Validated', {variant: 'success'}));
    } catch (e) {
      editor.diagnostics.value = BindingDialog.issues(e);
      this._validated.value = false;
      this._say(DomainErrors.message(e), true);
    }
    this._renderIssues();
  }

  /** The schema exists once the create answered: the promise resolves there, whatever the access
   * rows and the app opening afterwards report. */
  private async _create(): Promise<boolean> {
    const editor = this.editor!;
    const plan = editor.plan();
    this._say('Creating…');
    try {
      await grok.dapi.domains.createSchema(plan.name, {friendlyName: plan.friendlyName || undefined,
        manifest: plan.manifest});
    } catch (e) {
      editor.diagnostics.value = BindingDialog.issues(e);
      this._validated.value = false;
      this._renderIssues();
      this._say(DomainErrors.message(e), true);
      return false;
    }
    const problems = await this._applyAccess(plan);
    grok.dapi.domains.invalidateUiCaches();
    notify.info(`Domain schema ${plan.name} created`);
    for (const problem of problems)
      notify.warning(problem);
    this._settle(plan.name);
    try {
      const view = await route(`/domains/${plan.name}/${Object.keys(plan.manifest.tables)[0]}`);
      if (view !== null)
        grok.shell.addView(view);
    } catch (e) {
      notify.error(`Domain schema ${plan.name} could not be opened: ${DomainErrors.message(e)}`);
    }
    return true;
  }

  /** The access rows after the create: schema-scope grants fan out on the server, table grants
   * are row access, a visibility row restricts the column to its groups. Every failure is
   * collected and named — the schema exists by now, and the rows can be redone from its page. */
  private async _applyAccess(plan: ManifestPlan): Promise<string[]> {
    const problems: string[] = [];
    const groupId = (name: string): string => {
      const id = this._groups.get(name);
      if (id === undefined)
        throw new Error(`Group "${name}" is not known`);
      if (id === null)
        throw new Error(`Group "${name}" is ambiguous — two groups carry the name`);
      return id;
    };
    const domains = grok.dapi.domains;
    for (const g of plan.grants) {
      const permissions = [...(g.view ? ['View'] : []), ...(g.edit ? ['Edit'] : []), ...(g.delete ? ['Delete'] : [])];
      for (const permission of permissions) {
        try {
          const id = groupId(g.group);
          await (g.table === null ? domains.schema(plan.name).grant(id, permission) :
            domains.table(`${plan.name}.${g.table}`).grant(id, permission));
        } catch (e) {
          problems.push(`${g.table ?? plan.name}: ${permission} for ${g.group}: ${DomainErrors.message(e)}`);
        }
      }
    }
    for (const r of plan.restrictions) {
      const table = domains.table(`${plan.name}.${r.table}`);
      try {
        if (r.groups.length === 0)
          await table.restrictColumn(r.column);
        for (const group of r.groups)
          await table.shareColumn(r.column, groupId(group), 'View');
      } catch (e) {
        problems.push(`${r.table}.${r.column}: visibility: ${DomainErrors.message(e)}`);
      }
    }
    return problems;
  }
}

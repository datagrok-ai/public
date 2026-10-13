/* `domains.authoring.editBinding` — the "Edit binding" dialog over a REGISTERED external binding:
   **Connection** (the binding's connection, remote schema and catalog, read-only — a binding is
   never re-pointed), **Design** (the `ManifestEditor` in edit mode over the registered manifest,
   the access snapshot and a fresh draft of the warehouse), **Review** (the change list, each line
   annotated from the dry run — what a removal takes with it, what is already so; VALIDATE is the
   dry run bound to the exact apply body it checked, SAVE the one apply, access changes included)
   and **Saved** (what was applied — or that the server found nothing to change — OPEN). A `version-conflict` (an `access-conflict` too) reloads
   the binding and replays the edits over it, conflicts and drops reported for an explicit look on
   Design; a save whose answer was lost is confirmed through the registry, never replayed — one
   the registry cannot confirm ends the edit, as does a binding deleted and re-created under the
   name: REOPEN starts over from the server's version. The server owns every rule — the dialog
   names its refusals. */
import * as grok from 'datagrok-api/grok';
import type * as DG from 'datagrok-api/dg';
import {computed, signal} from '../../../core/signals.js';
import {div, divV, span} from '../../../core/elements.js';
import {plural} from '../../../core/text.js';
import {Section} from '../../../components/containers/section.js';
import {BoolInput} from '../../../components/inputs/bool-input.js';
import {badge} from '../../../components/display/badge.js';
import {DomainErrors} from '../errors.js';
import {groupInput} from '../../inputs/group-input.js';
import {BindingWizard} from './binding-dialog.js';
import {ManifestEditor} from './manifest-editor.js';
import type {ApplyPayload} from './manifest-editor.js';
import {ManifestModel} from './manifest-model.js';
import type {AccessPrincipal, AccessSnapshot, DraftEnvelope, ManifestChange, ManifestJson, ManifestTableJson,
  RebaseOp, RebaseReport} from './manifest-model.js';

/** The platform's Dart client learns of an apply from the apply itself; one confirmed through the
 * registry after a lost answer is announced to it here (read at call time like the create's). */
const api = globalThis as {grok_Dapi_Domains_SchemaAltered?: (dart: unknown, name: string) => void};

export interface EditBindingResult {
  name: string;
  /** The apply's answer — the plan with the access effects as resolved under the lock, `applied:
   * false` with `noop: true` where the server found nothing to change; null where the answer was
   * lost and the registry confirmed the save. */
  applied: DG.DomainApplied | null;
}

/** Opens the dialog; resolves once it closes — to the saved binding with the apply's answer, null
 * where it was cancelled before the save. Refuses while domain databases are off, and — by the
 * registry's own refusal — over a schema the caller may not read: the dialog does not open. */
export async function editBinding(schema: string): Promise<EditBindingResult | null> {
  if (grok.shell.settings.enableDomainDatabases !== true)
    throw new Error('Domain databases are a Beta feature — enable them in Settings > Beta');
  return new BindingEditDialog(schema).open();
}

/** Everything an edit reads about the binding — at open, and again after a conflict. */
interface Loaded {
  manifest: ManifestJson;
  snapshot: AccessSnapshot | null;
  draft: DraftEnvelope | null;
  friendlyName: string;
  description: string;
  /** The reads that failed, by name — the status line says so; the model shows them as unknown. */
  problems: string[];
}

/** The registry read the dialog waits for before it shows, and the reads that go on behind it. */
interface Read {
  manifest: ManifestJson;
  rest: Promise<Omit<Loaded, 'manifest'>>;
}

export class BindingEditDialog extends BindingWizard<EditBindingResult> {
  private readonly _name: string;
  private readonly _handle: DG.DomainSchemaClient;
  /** Why the editor could not be built. */
  private readonly _problem = signal<string | null>(null);
  /** The dry run's plan of the validated payload — what annotates the change list. */
  private readonly _plan = signal<DG.DomainApplyPlan | null>(null);
  /** Why the edit ended short of a save — its outcome unknown, or the binding another schema
   * now — with the facts; REOPEN is the one way on. */
  private readonly _ended = signal<{why: string, facts: string[]} | null>(null);
  private readonly _facts = divV([], 'u2-binding-connection u2-binding-edit-facts');
  private readonly _designHost = div([], 'u2-binding-design');
  private readonly _changes = divV([], 'u2-binding-changes');
  private readonly _body = document.createElement('pre');
  private readonly _outcome = divV([], 'u2-binding-reloaded u2-binding-outcome');
  private readonly _reloaded = divV([], 'u2-binding-reloaded');
  private readonly _confirm: BoolInput;
  private readonly _applied = divV([], 'u2-binding-changes');
  private readonly _savedHost = divV([], 'u2-binding-created u2-binding-saved');
  /** False after a reload that conflicted or dropped an edit, until Design was looked at again. */
  private _looked = true;
  private _report: RebaseReport | null = null;
  private _result: EditBindingResult | undefined;
  /** The dialog REOPEN started over: what the caller gets is its outcome. */
  private _reopened: Promise<EditBindingResult | null> | undefined;
  private _saved: ManifestChange[] = [];
  private _savedVersion = '';
  private _firstTable = '';

  constructor(schema: string) {
    super();
    this._name = schema;
    this._handle = grok.dapi.domains.schema(schema);
    this.root.dataset.u2 = 'binding-edit-dialog';
    this._body.className = 'u2-binding-json';
    this._confirm = this.runInScope(() => new BoolInput({label: 'I understand what is removed', name: 'confirmDestructive'}));
    this._confirm.root.classList.add('u2-binding-confirm');
    this._confirm.root.style.display = 'none';
    this._mount({
      start: 'design',
      steps: [
        {id: 'connection', title: 'Connection', content: this._facts},
        {id: 'design', title: 'Design', content: () => this._designHost, canProceed: () => this._designGate(),
          onActivate: () => this._looked = true},
        {id: 'review', title: 'Review', content: () => this._reviewStep(),
          nextText: computed(() => this._ended.value === null ? 'SAVE' : 'REOPEN'),
          onActivate: () => this._review(), canProceed: () => this._reviewGate(),
          commit: () => this._ended.peek() === null ? this._save() : this._reopen(),
          actions: [{text: 'VALIDATE', run: () => this._validate(), visible: computed(() => this._ended.value === null),
            enabled: computed(() => this._ended.value === null && this._editor.value !== undefined &&
              this._editor.value.model.blockers.value.length === 0)}]},
        {id: 'saved', title: 'Saved', content: () => this._savedHost, done: true,
          actions: [{text: 'OPEN', run: () => this._openAndClose()}]},
      ],
      onFinish: () => this._settle(this._result ?? null),
      onCancel: () => this._settle(null),
    });
  }

  /** The registry is read before the dialog shows: a schema it refuses is a refusal here, by
   * name. The other reads go on behind Design, which waits on them. A dialog REOPEN started
   * over answers for this one. */
  async open(): Promise<EditBindingResult | null> {
    return this._show(await this._read());
  }

  private async _read(): Promise<Read> {
    const {manifest, rest} = this._start();
    try {
      return {manifest: await manifest, rest};
    } catch (e) {
      this.dispose();
      throw new Error(`Binding ${this._name} could not be read: ${BindingEditDialog._failure(e)}`);
    }
  }

  private _show(read: Read): Promise<EditBindingResult | null> {
    const settled = this._open(`Edit binding ${this._name}`);
    void this._build(read.rest.then((rest) => ({manifest: read.manifest, ...rest})));
    return settled.then((result) => this._reopened ?? result);
  }

  /** The apply body as one string, keys sorted — what a validation is bound to; the access deltas
   * are part of it. */
  protected _payloadOf(editor: ManifestEditor): string {
    const model = editor.model;
    model.revision.value;
    model.friendlyName.value;
    model.description.value;
    editor.access.grants.value;
    editor.access.visibility.value;
    return BindingEditDialog._canonical(editor.editPlan().payload);
  }

  /** The reads: the manifest is what the dialog waits for, the rest settle into what the model
   * can show as unknown — each failure by name. The warehouse is asked only once the registry
   * answered; the registry reads go at once. */
  private _start(): {manifest: Promise<ManifestJson>, rest: Promise<Omit<Loaded, 'manifest'>>} {
    const manifest = this._handle.manifest() as Promise<ManifestJson>;
    const snapshot = BindingEditDialog._within(this._handle.access(), 'The access snapshot');
    const draft = manifest.then(() => BindingEditDialog._within(this._handle.draft(), 'The warehouse catalog'), () => null);
    const entity = BindingEditDialog._entity(this._name);
    const problems: string[] = [];
    const settle = async <T>(read: Promise<T>, what: string): Promise<T | null> => {
      try {
        return await read;
      } catch (e) {
        problems.push(`${what} could not be read: ${BindingEditDialog._failure(e)}`);
        return null;
      }
    };
    const rest = (async (): Promise<Omit<Loaded, 'manifest'>> => {
      const [s, d, e] = await Promise.all([settle(snapshot, 'The access snapshot'),
        settle(draft, 'The warehouse catalog'), settle(entity, 'The schema\'s caption')]);
      return {snapshot: s as AccessSnapshot | null, draft: d as DraftEnvelope | null,
        friendlyName: e?.friendlyName ?? '', description: e?.description ?? '', problems};
    })();
    return {manifest, rest};
  }

  private async _reload(): Promise<Loaded> {
    const {manifest, rest} = this._start();
    const m = await manifest;
    return {manifest: m, ...await rest};
  }

  /** The editor over what was read; a snapshot and a manifest read at different versions are
   * read once more, and said if they still disagree. */
  private async _build(loading: Promise<Loaded>): Promise<void> {
    this._say('Reading the binding…');
    try {
      let loaded = await loading;
      let editor = this._editorOver(loaded);
      if (editor.stale) {
        editor.dispose();
        loaded = await this._reload();
        editor = this._editorOver(loaded);
      }
      if (this.scope.isDisposed) {
        editor.dispose();
        return;
      }
      this._renderFacts(loaded.manifest);
      this._diagnosed = this._payloadOf(editor);
      this._designHost.replaceChildren(editor.root);
      this._editor.value = editor;
      this._say(BindingEditDialog._notes(loaded, editor).join('; '), loaded.problems.length > 0);
    } catch (e) {
      this._problem.value = BindingEditDialog._failure(e);
      this._say(this._problem.value, true);
    }
  }

  private _editorOver(loaded: Loaded): ManifestEditor {
    return this.runInScope(() => new ManifestEditor(loaded.draft, {
      context: {mode: 'edit', storage: 'external'}, baseline: loaded.manifest, snapshot: loaded.snapshot,
      friendlyName: loaded.friendlyName, description: loaded.description, author: BindingEditDialog._author(),
      principalPicker: (onPick) => groupInput({onPick: (g, label) => onPick({id: g.id, label})}).root,
    }));
  }

  /** The author's own personal group: restricting a column keeps them on it, explicitly. */
  private static _author(): AccessPrincipal | undefined {
    const user = grok.shell.user;
    const group = user?.group;
    return group == null ? undefined : {id: group.id, label: user.friendlyName};
  }

  private static _notes(loaded: Loaded, editor: ManifestEditor): string[] {
    const notes = [...loaded.problems];
    if (editor.stale)
      notes.push('The access snapshot and the manifest were read at different versions — the access rows may be off');
    return notes;
  }

  private _renderFacts(manifest: ManifestJson): void {
    const storage = manifest.storage;
    const fact = (label: string, value: unknown): HTMLElement =>
      div([span(label, 'u2-binding-fact-label'), span(value === undefined || value === null || value === '' ? '—' : String(value))],
        'u2-binding-fact');
    this._facts.replaceChildren(fact('Connection', storage?.connection), fact('Remote schema', storage?.schema),
      fact('Catalog', storage?.catalog),
      span('A binding is never re-pointed: register a new one for another connection or schema.', 'u2-binding-facts'));
  }

  private _designGate(): string | null {
    const problem = this._problem.value;
    if (problem !== null)
      return problem;
    const editor = this._editor.value;
    if (editor === undefined || this._current.value === null)
      return 'Reading the binding…';
    const plan = editor.editPlan();
    if (plan.blockers.length > 0) {
      return `${plural(plan.blockers.length, 'registered item', 'registered items')} cannot bind as declared — ` +
        'take them out, or fix the warehouse';
    }
    return plan.changes.length === 0 ? 'Nothing changed' : null;
  }

  private _reviewStep(): HTMLElement {
    const body = this.runInScope(() => new Section({title: 'Apply body', expanded: false}));
    body.root.classList.add('u2-binding-body');
    body.body.append(this._body);
    return div([divV([this._changes, body.root], 'u2-binding-changes-pane'),
      divV([this._outcome, this._reloaded, this._confirm.root, this._issues], 'u2-binding-review-side')], 'u2-binding-review');
  }

  private _reviewGate(): string | null {
    if (this._ended.value !== null)
      return null;
    if (!this._looked)
      return 'Look over the reloaded edits on Design';
    if (!this._isValidated())
      return 'Validate before saving';
    return this._plan.value?.destructive === true && !this._confirm.value.value ? 'Confirm what is removed' : null;
  }

  /** Rebuilt on every activation and after every validation: the change list annotated from the
   * plan while the validation holds, the findings, the reload report, the confirmation — or,
   * once the edit ended, why. */
  private _review(): void {
    const editor = this.editor!;
    const {payload, changes} = editor.editPlan();
    const validated = this._isValidated();
    const plan = validated ? this._plan.peek() : null;
    this._body.textContent = JSON.stringify(payload, null, 2);
    const ended = this._ended.peek();
    this._renderChanges(this._changes, changes, plan,
      ended !== null ? 'Changes' : validated ? 'Changes, as validated' : 'Changes, not validated');
    this._confirm.root.style.display = plan?.destructive === true ? '' : 'none';
    this._renderOutcome();
    this._renderReloaded();
    this._say(ended !== null ? `${ended.why} — see Review` : validated ? badge('Validated', {variant: 'success'}) : '', ended !== null);
    this._renderIssues();
    this._issues.style.display = ended !== null ? 'none' : '';
  }

  private _renderOutcome(): void {
    const ended = this._ended.peek();
    if (ended === null) {
      this._outcome.replaceChildren();
      return;
    }
    this._outcome.replaceChildren(span(ended.why, 'u2-binding-changes-title'),
      ...ended.facts.map((fact) => div([span(fact)], 'u2-binding-change')),
      span('Nothing was replayed. REOPEN starts over from the server\'s version; the edits above are not carried over.',
        'u2-binding-change-note'));
  }

  /** Ends the edit: no save, no rebase; REOPEN in place of SAVE, and Design says so. */
  private _end(why: string, facts: string[]): void {
    this._ended.value = {why, facts};
    this._setPlan(null);
    this._validated.value = null;
    this._designHost.prepend(div([span(`${why}: nothing here will be saved — REOPEN on Review starts over from the server's version`)],
      'u2-manifest-panel-note u2-binding-ended'));
    this._review();
  }

  /** The dialog over the server's version takes this one's place once it has read the registry;
   * a read that fails leaves this one standing, and says why. */
  private async _reopen(): Promise<false> {
    this._say('Reading the binding again…');
    const next = new BindingEditDialog(this._name);
    let read: Read;
    try {
      read = await next._read();
    } catch (e) {
      this._say(`${(e as Error).message} — REOPEN again`, true);
      return false;
    }
    this._reopened = next._show(read);
    this._settle(null);
    return false;
  }

  private _renderChanges(host: HTMLElement, changes: ManifestChange[], plan: DG.DomainApplyPlan | null, title: string): void {
    host.replaceChildren(span(title, 'u2-binding-changes-title'));
    for (const change of changes) {
      const note = plan === null ? '' : BindingEditDialog._annotate(change, plan);
      host.append(div([span(change.text), ...(note === '' ? [] : [span(note, 'u2-binding-change-note')])],
        `u2-binding-change${change.removes ? ' u2-binding-change-removes' : ''}`));
    }
    for (const note of plan === null ? [] : BindingEditDialog._planNotes(plan))
      host.append(span(note, 'u2-binding-change-note'));
  }

  /** What the plan says about one line of the change list, by its id: what a removal takes with
   * it (`lost`), what a renamed caption was (`metadata`), an access op that is already so. */
  private static _annotate(change: ManifestChange, plan: DG.DomainApplyPlan): string {
    const [kind, ...rest] = change.id.split(':');
    if (kind === 'table' && change.removes) {
      const lost = plan.lost?.tables[change.table!];
      return lost === undefined ? '' : BindingEditDialog._lostTable(lost);
    }
    if (kind === 'column' && change.removes) {
      const lost = plan.lost?.columns[`${change.table}.${change.column}`];
      return lost === undefined ? '' : BindingEditDialog._lostColumn(lost);
    }
    if (kind === 'schema') {
      const m = plan.metadata?.[rest[0]];
      return m === undefined ? '' : `was ${m.from === null || m.from === '' ? 'empty' : `"${m.from}"`}`;
    }
    return kind === 'access' && plan.access !== undefined ? BindingEditDialog._effect(rest, plan.access) : '';
  }

  private static _lostTable(lost: DG.DomainLostTable): string {
    const parts: string[] = [];
    const grants = lost.grants + lost.coreSchemaGrants;
    const restricted = Object.keys(lost.restrictions).length;
    if (grants > 0)
      parts.push(plural(grants, 'grant', 'grants'));
    if (restricted > 0)
      parts.push(plural(restricted, 'restricted column', 'restricted columns'));
    if (lost.promotedRows > 0)
      parts.push(`${plural(lost.promotedRows, 'promoted row', 'promoted rows')} with ${plural(lost.rowGrants, 'grant', 'grants')}`);
    if (lost.savedFilters > 0)
      parts.push(plural(lost.savedFilters, 'saved filter', 'saved filters'));
    const one = parts.length === 1 && grants + restricted + lost.promotedRows + lost.savedFilters === 1;
    let text = parts.length === 0 ? 'nothing else goes with it' : `${parts.join(', ')} ${one ? 'goes' : 'go'} with it`;
    if (lost.affectedFilters > 0)
      text += `; ${BindingEditDialog._unresolved(lost.affectedFilters, 'of other tables')}`;
    return `${text}; the warehouse is untouched`;
  }

  private static _lostColumn(lost: DG.DomainLostColumn): string {
    let text = !lost.restricted ? 'nothing else goes with it' : lost.grants === 0 ? 'its restriction goes with it' :
      lost.grants === 1 ? 'its restriction and the 1 grant on it go with it, whoever holds it' :
        `its restriction and all ${lost.grants} grants on it go with it, whoever holds them`;
    if (lost.affectedFilters > 0)
      text += `; ${BindingEditDialog._unresolved(lost.affectedFilters, 'naming it')}`;
    return `${text}; the warehouse is untouched`;
  }

  private static _unresolved(filters: number, which: string): string {
    return `${plural(filters, 'saved filter', 'saved filters')} ${which} no longer ${filters === 1 ? 'resolves' : 'resolve'}`;
  }

  /** The effect of the access op behind `access:<op>:<target>…`: "already so" where the server
   * found the triple there (or gone) already. */
  private static _effect(id: string[], access: DG.DomainAccessEffects): string {
    const [op, target, group, permission] = id;
    if (op === 'grant' || op === 'revoke') {
      const effect = access[op].find((e) => e.table === target && e.group.id === group && e.permission === permission);
      return effect?.effect === 'none' ? 'already so' : '';
    }
    const dot = target.indexOf('.');
    const table = target.slice(0, dot);
    const column = target.slice(dot + 1);
    if (op === 'unrestrict') {
      const effect = access.unrestrict.find((e) => e.table === table && e.column === column);
      return effect?.effect === 'none' ? 'already visible to everyone' : '';
    }
    const effect = access.restrict.find((e) => e.table === table && e.column === column);
    if (effect === undefined)
      return '';
    const parts: string[] = [];
    if (effect.effect === 'none')
      parts.push('already restricted');
    for (const g of effect.grants.filter((x) => x.effect === 'none'))
      parts.push(`${g.permission} for ${g.group.friendlyName}: already so`);
    for (const g of effect.revokes.filter((x) => x.effect === 'none'))
      parts.push(`${g.permission} revoked from ${g.group.friendlyName}: already so`);
    return parts.join('; ');
  }

  /** What the plan says beyond the change list: key columns the apply makes visible again, and
   * refs whose saved filters stop resolving. */
  private static _planNotes(plan: DG.DomainApplyPlan): string[] {
    const notes: string[] = [];
    if (plan.keyColumnsUnrestricted.length > 0) {
      notes.push(`Key columns visible to everyone again: ${plan.keyColumnsUnrestricted
        .map((c) => `${c.table}.${c.column}`).join(', ')}`);
    }
    for (const [address, ref] of Object.entries(plan.lost?.refs ?? {}))
      notes.push(`${address}: ${BindingEditDialog._unresolved(ref.affectedFilters, 'through it')}`);
    return notes;
  }

  private _renderReloaded(): void {
    const report = this._report;
    if (report === null) {
      this._reloaded.replaceChildren();
      return;
    }
    const line = (op: RebaseOp, tail: string): HTMLElement => div([span(op.text), span(tail, 'u2-binding-change-note')],
      'u2-binding-change');
    const value = (v: unknown): string => v === null || v === undefined || v === '' ? 'empty' : JSON.stringify(v);
    this._reloaded.replaceChildren(
      span(`Reloaded: ${plural(report.applied.length, 'edit', 'edits')} kept, ` +
        `${plural(report.conflicts.length, 'conflict', 'conflicts')}, ${report.dropped.length} dropped`, 'u2-binding-changes-title'),
      ...report.conflicts.map((op) => line(op, op.note ?? `the server's ${value(op.server)} stands; yours was ${value(op.to)}`)),
      ...report.dropped.map((op) => line(op, 'no longer there')));
  }

  /** The dry run over the apply body as it stands; an answer to a payload since changed, or to an
   * older run, is dropped. A failure the server did not answer is no finding. A destructive plan
   * the server hands back with its confirmation demand is the plan all the same. */
  private async _validate(): Promise<void> {
    const editor = this.editor!;
    const {payload} = editor.editPlan();
    const key = this._payloadOf(editor);
    const gen = ++this._validateGen;
    this._say('Validating…');
    let refusal: unknown = null;
    let plan: DG.DomainApplyPlan | null = null;
    try {
      plan = await this._handle.apply(payload, {dryRun: true});
    } catch (e) {
      const demanded = (e as {body?: {plan?: DG.DomainApplyPlan}} | null)?.body?.plan;
      if (DomainErrors.codeOf(e) === 'destructive-confirmation-required' && demanded !== undefined)
        plan = demanded;
      else
        refusal = e;
    }
    if (gen !== this._validateGen)
      return;
    if (this.editor !== editor || this._payloadOf(editor) !== key)
      return this._say('Changed while validating — validate again');
    if (refusal !== null && !BindingEditDialog._answered(refusal))
      return this._say(BindingEditDialog._failure(refusal), true);
    this._diagnosed = key;
    editor.diagnostics.value = refusal === null ? [] : BindingEditDialog.issues(refusal);
    this._validated.value = refusal === null ? key : null;
    this._setPlan(plan);
    this._review();
    if (refusal !== null)
      this._say(DomainErrors.message(refusal), true);
  }

  /** A confirmation holds for the plan it was given: another plan wants its own. */
  private _setPlan(plan: DG.DomainApplyPlan | null): void {
    this._plan.value = plan;
    this._confirm.value.value = false;
  }

  /** The one apply. A conflict reloads the binding and replays the edits; a lost answer asks the
   * registry whether the save landed; any other refusal lands on the findings, SAVE still offered. */
  private async _save(): Promise<boolean> {
    const editor = this.editor!;
    const {payload, changes} = editor.editPlan();
    const body = this._plan.peek()?.destructive === true ? {...payload, confirmDestructive: true} : payload;
    this._say('Saving…');
    let answer: DG.DomainApplied | null;
    try {
      answer = await this._handle.apply(body);
    } catch (e) {
      const code = DomainErrors.codeOf(e);
      if (code === 'version-conflict' || code === 'access-conflict') {
        await this._rebase(payload, e);
        return false;
      }
      if (BindingEditDialog._answered(e)) {
        this._diagnosed = this._payloadOf(editor);
        editor.diagnostics.value = BindingEditDialog.issues(e);
        this._validated.value = null;
        this._renderIssues();
        this._say(DomainErrors.message(e), true);
        return false;
      }
      if (!await this._landed(payload, e))
        return false;
      answer = null;
      api.grok_Dapi_Domains_SchemaAltered?.(grok.dapi.domains.dart, this._name);
    }
    this._result = {name: this._name, applied: answer};
    this._saved = changes;
    this._firstTable = Object.keys(editor.model.toJSON().tables)[0] ?? '';
    grok.dapi.domains.invalidateUiCaches();
    this._renderSaved();
    return true;
  }

  /** Whether an apply that got no answer landed cleanly: the registry is one version past the
   * edited one, in the same incarnation, AND holds this edit. Anything else — the version moved
   * on, the next version is someone else's, the name is another schema's, the reads fail — is an
   * outcome nobody saw, and ends the edit rather than replay it. */
  private async _landed(payload: ApplyPayload, e: unknown): Promise<boolean> {
    this._say('No answer to the save — reading the registry…');
    const next = String(Number(payload.ifVersion) + 1);
    let fact: string;
    try {
      const manifest = await this._handle.manifest();
      const same = payload.ifIncarnation === undefined || manifest.incarnation === payload.ifIncarnation;
      if (same && manifest.version === next && await this._holds(payload, manifest)) {
        this._savedVersion = next;
        return true;
      }
      fact = !same ? BindingEditDialog._recreated(this._name, manifest.version) :
        manifest.version === next ? `version ${next} of ${this._name} holds someone else's save, not this edit` :
          `${this._name} is at version ${manifest.version} (this edit was made against ${payload.ifVersion})`;
    } catch (x) {
      fact = `the registry could not be read (${BindingEditDialog._failure(x)})`;
    }
    // it may have landed: whoever shows the binding reads it again
    grok.dapi.domains.invalidateUiCaches();
    api.grok_Dapi_Domains_SchemaAltered?.(grok.dapi.domains.dart, this._name);
    this._end('Outcome unknown', [`No answer to the save (${BindingEditDialog._failure(e)}) — whether it landed is unknown`, fact]);
    return false;
  }

  /** Whether the registry holds this edit: every table sent as registered, none dropped still
   * there, the metadata and writability as sent, and every access op's effect in place — read
   * at the manifest's version; a grant the snapshot cannot show proves nothing. */
  private async _holds(payload: ApplyPayload, manifest: DG.DomainRegisteredManifest): Promise<boolean> {
    const canonical = (decl: ManifestTableJson, logical: string): string =>
      BindingEditDialog._canonical(ManifestModel.canonicalTable(decl, logical));
    for (const [logical, decl] of Object.entries(payload.tables ?? {})) {
      const registered = manifest.tables[logical] as ManifestTableJson | undefined;
      if (registered === undefined || canonical(registered, logical) !== canonical(decl, logical))
        return false;
    }
    if ((payload.dropTables ?? []).some((t) => manifest.tables[t] !== undefined))
      return false;
    if (payload.storage !== undefined && (manifest.storage?.writable ?? false) !== payload.storage.writable)
      return false;
    if (payload.friendlyName !== undefined || payload.description !== undefined) {
      const entity = await BindingEditDialog._entity(this._name);
      if (entity === null || (payload.friendlyName !== undefined && entity.friendlyName !== payload.friendlyName) ||
          (payload.description !== undefined && (entity.description ?? '') !== payload.description))
        return false;
    }
    const access = payload.access;
    if (access === undefined)
      return true;
    const snapshot = await this._handle.access();
    if (snapshot.version !== manifest.version || snapshot.incarnation !== manifest.incarnation)
      return false;
    const held = (table: string, group: string, permission: string): boolean | null => {
      const grants = snapshot.tables[table]?.grants;
      return grants == null ? null : grants.some((g) => g.group.id === group && g.permissions.includes(permission as DG.DomainPermission));
    };
    if (access.grant.some((op) => held(op.table, op.group, op.permission) !== true) ||
        access.revoke.some((op) => held(op.table, op.group, op.permission) !== false))
      return false;
    for (const op of access.restrict) {
      const column = snapshot.columns[`${op.table}.${op.column}`];
      if (column?.state !== 'restricted')
        return false;
      const on = (t: {group: string, permission: string}): boolean | null => {
        const groups = t.permission === 'Edit' ? column.edit : column.view;
        return groups == null ? null : groups.some((g) => g.id === t.group);
      };
      if (op.grant.some((t) => on(t) !== true) || op.revoke.some((t) => on(t) !== false))
        return false;
    }
    return access.unrestrict.every((op) => snapshot.columns[`${op.table}.${op.column}`]?.state === 'unrestricted');
  }

  private static _recreated(name: string, version: string | undefined): string {
    return `${name} was deleted and re-created since — another schema is registered under the name now, at version ${version ?? '?'}`;
  }

  /** After a `version-conflict`: the binding as it is now, the edits replayed over it. What
   * conflicted or vanished is reported on Review and wants a look at Design before a save. A
   * name that is another schema's now (deleted and re-created since) takes no replay. */
  private async _rebase(payload: ApplyPayload, e: unknown): Promise<void> {
    const editor = this.editor!;
    const conflict = e as {currentIncarnation?: unknown, currentVersion?: unknown};
    if (typeof conflict.currentIncarnation === 'string' && payload.ifIncarnation !== undefined &&
        conflict.currentIncarnation !== payload.ifIncarnation)
      return this._end('The binding you edited is gone', [BindingEditDialog._recreated(this._name, String(conflict.currentVersion ?? '?'))]);
    this._say('The binding changed meanwhile — reloading…');
    let loaded: Loaded;
    try {
      loaded = await this._reload();
    } catch (x) {
      this._say(`The binding changed meanwhile, and could not be re-read: ${BindingEditDialog._failure(x)} — SAVE again`, true);
      return;
    }
    if (loaded.manifest.incarnation !== editor.model.incarnation)
      return this._end('The binding you edited is gone', [BindingEditDialog._recreated(this._name, loaded.manifest.version)]);
    const report = editor.rebase(loaded.manifest, loaded.snapshot, loaded.draft,
      {friendlyName: loaded.friendlyName, description: loaded.description});
    this._report = report;
    this._setPlan(null);
    this._validated.value = null;
    this._diagnosed = this._payloadOf(editor);
    this._looked = report.conflicts.length === 0 && report.dropped.length === 0;
    this._review();
    const notes = [`Reloaded at version ${loaded.manifest.version ?? '?'}: ${plural(report.applied.length, 'edit', 'edits')} kept, ` +
      `${plural(report.conflicts.length, 'conflict', 'conflicts')}, ${report.dropped.length} dropped — validate again`,
    ...BindingEditDialog._notes(loaded, editor)];
    this._say(notes.join('; '), !this._looked || loaded.problems.length > 0);
  }

  private _renderSaved(): void {
    const answer = this._result!.applied;
    if (answer?.noop === true) {
      const nothing = 'Nothing to save — the binding already holds this';
      this._savedHost.replaceChildren(div([badge('Unchanged'),
        span(`${nothing}; domain schema ${this._name} stays at version ${answer.version}`)], 'u2-binding-created-head'));
      this._say(nothing);
      return;
    }
    const version = answer?.version ?? this._savedVersion;
    const head = div([badge('Saved', {variant: 'success'}),
      span(answer === null ? `Domain schema ${this._name} is at version ${version} — the answer was lost, the registry confirms it` :
        `Domain schema ${this._name} is at version ${version}`)], 'u2-binding-created-head');
    this._savedHost.replaceChildren(head, this._applied);
    this._renderChanges(this._applied, this._saved, answer, 'Applied');
    this._say(badge('Saved', {variant: 'success'}));
  }

  private async _openAndClose(): Promise<void> {
    const result = this._result!;
    const first = this._firstTable;
    this._settle(result);
    await BindingEditDialog._openApp(result.name, first);
  }
}

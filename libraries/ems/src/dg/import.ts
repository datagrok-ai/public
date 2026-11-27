import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import type * as DG from 'datagrok-api/dg';
import {signal, computed, Control, div, divH, divV, span, plural, text, link, Wizard, DataTable, ChoiceInput,
  BoolInput, badge, dot, statusLine, Form, notify, Scope, isEmpty} from '@datagrok-libraries/u2';
import type {IProperty, ChoiceItem, BadgeVariant, RowView} from '@datagrok-libraries/u2';
import {capital, nounOf} from '../core/domain-backend.js';
import {MemoryEditState} from '../core/edit-state.js';
import type {DomainBatchOptionsLike, DomainBatchReportLike, DomainBatchRowLike, DomainBatchValidationLike}
  from '../core/domain-backend.js';
import {fromDartInput} from '@datagrok-libraries/u2/dg';
import {DomainTable} from './facade.js';
import {DomainDataTable} from './grid.js';
import {DomainErrors} from './errors.js';
import {ImportLinks} from './import-links.js';
import {DomainAddress} from './address.js';
import type {ImportLink, LinkField, LinkMatch} from './import-links.js';

type Placed = DomainBatchRowLike & {index: number, id: string | null};
type Report = Omit<DomainBatchReportLike, 'rows'> & {rows: Placed[]};

export interface DomainImportOptions {
  /** The frame to import, when the caller already has one; the picker still offers the rest. */
  source?: DG.DataFrame;
}

/** Everything an app can make its own about importing its table — `table.importer`. Every member
 * is optional: without any, the import maps by name and semantic type, previews the server's dry
 * run, and reports every row's outcome and key. */
export interface DomainImporter {
  /** The columns {@link prepare} fills: the mapping does not offer them. */
  derives?: string[];
  /** The target a source column goes to, asked before the name and semantic-type matching;
   * undefined leaves the column to it. */
  match?(column: DG.Column, targets: IProperty[]): string | undefined;
  /** One verdict per row (the mapped cells under their target names), in order; `progress` takes
   * the number of rows done. */
  prepare?(rows: Record<string, unknown>[], progress: (done: number) => void): Promise<DomainImportVerdict[]>;
  /** The columns of the preview or the report, given the ones the import would show: `~row`,
   * `~result`, `~note` and each mapped target by its name in the preview; `~row`, `~result`, `~key` and
   * `~note` in the report. */
  columns?(defaults: DomainImportColumn[], step: DomainImportStep): DomainImportColumn[];
  /** The height of a preview or report row, for columns that draw. */
  rowHeight?: number;
  /** How a refusal reads in a row's note; null leaves it to the default ("Title: Value can't be
   * empty", "Code is already used"). */
  describe?(error: DomainImportError, caption: string): string | null;
  /** The source column each row's key is written back to, offered as an option; false offers none.
   * By default the business key's caption, where the key is one column. */
  keyColumn?: string | false;
}

export type DomainImportStep = 'preview' | 'report';

/** What becomes of one row: it is sent with `values` added, or kept out as the row the table
 * already has (`existing`, its id), as a repeat of an earlier row of the same file (`sameAs`, its
 * index) or as an `error`. `label` names it on the row's badge where it has no verb, `note` says why. */
export interface DomainImportVerdict {
  values?: Record<string, unknown>;
  existing?: string;
  sameAs?: number;
  error?: string;
  label?: string;
  note?: string;
}

/** A refusal of one row, by the column it is about where it names one. */
export interface DomainImportError {
  column?: string;
  code?: string;
  message: string;
}

/** A source row as the preview and the report show it. */
export interface DomainImportRow {
  /** Its place in the source frame, from 0. */
  index: number;
  /** The mapped cells under their target names, with what the importer derived. */
  values: Record<string, unknown>;
  /** What the batch will do or did with it — `insert`, `update`, `skip`, `error`; `aborted` where an
   * all-or-nothing import committed nothing; '' where nothing says. */
  kind: string;
  /** The badge text the verdict chose; '' leaves it to the kind. */
  label: string;
  note: string;
  errors: DomainImportError[];
  /** The row it is now — in the report, once read back — by id, and by its key: the key column's value where the
   * business key is one column, else the row's caption. */
  id: string | null;
  key: string;
}

/** A column of the preview or the report. */
export interface DomainImportColumn {
  name: string;
  header: string;
  width?: string;
  render(row: DomainImportRow): HTMLElement | string;
}

/** A source row as the import plans it: the mapped cells, what it sends (null where the verdict
 * keeps it out), and the verdict. */
interface Planned {
  cells: Record<string, unknown>;
  payload: Record<string, unknown> | null;
  verdict: DomainImportVerdict;
}

/** What a source column maps to when it maps to nothing. */
const SKIP = '(skip)';
const SKIP_LABEL = 'Don\'t import';
/** The preview scans at most this many rows — it is advisory, and the server re-validates every
 * row it is sent; the commit itself is unbounded. */
const PREVIEW_ROWS = 1000;
/** Rows the preview draws — what the import looks like, not all of it. */
const PREVIEW_SHOWN = 20;
/** A row's kind as the preview shows it: the verb, not the wire's op name. */
const VERDICTS: Record<string, {label: string, variant: BadgeVariant}> = {
  insert: {label: 'Add', variant: 'success'},
  update: {label: 'Update', variant: 'success'},
  skip: {label: 'Skip', variant: 'default'},
  error: {label: 'Error', variant: 'error'},
};
/** A row's kind as the report shows it, once the batch has run. */
const OUTCOMES: Record<string, {label: string, variant: BadgeVariant}> = {
  insert: {label: 'Added', variant: 'success'},
  update: {label: 'Updated', variant: 'success'},
  skip: {label: 'Skipped', variant: 'default'},
  error: {label: 'Error', variant: 'error'},
  aborted: {label: 'Not imported', variant: 'default'},
};
/** Rows of the report's key lookup per query. */
const KEY_CHUNK = 500;
const USER_FILTER_BYTES = 6000;
const USER_READS = 4;

/** One row of the mapping step — the target is the choice input the flow keeps for that column. */
interface Mapping {column: string}

/** Opens the wizard over `table`; resolves to the server's report, or to null where it was
 * cancelled. Exported as `domains.import` — `import` cannot name a function. Only the mapped columns are posted,
 * under their target names, plus what the table's {@link DomainImporter} derives; the preview is the server's
 * dry run, the commit's own checks inside a rolled-back transaction. */
export async function openImport(table: DomainTable, options: DomainImportOptions = {}):
  Promise<DomainBatchReportLike | null> {
  if (table.table.batch === undefined) {
    notify.error(`${table.address}: the backend does not support batch import`);
    return null;
  }
  if (!table.table.support.transaction) {
    notify.error(`${table.address}: the table does not accept writes`);
    return null;
  }
  // an import inserts: the draft policy, so a key typed once on insert is a target too
  const writable = table.access.columnPolicy().forDraft();
  const derived = table.importer?.derives ?? [];
  const stuck = DomainTable.unsettable(table.properties, table.info, table.access, derived);
  const inserts = stuck.length === 0 && !table.info.hiddenRequired;
  const key = table.info.businessKey;
  const updates = table.table.support.batch.upsert && key.length > 0 &&
    key.every((k) => table.properties.some((p) => p.name === k)) && !stuck.some((p) => key.includes(p.name!));
  const auto = updates ? await autoNumbered(table) : null;
  const targets = table.properties.filter((p) => !derived.includes(p.name!) &&
    (writable.field(p.name!) === 'editable' || p.name === auto));
  if (targets.length === 0) {
    notify.warning(`There is no column of ${table.address} you may write.`);
    return null;
  }
  if (!inserts && !updates) {
    notify.warning(DomainTable.noRequired(table.info, stuck, table.access));
    return null;
  }
  // the wizard opens after the linked tables load: past that await it is a tree of its own
  const links = await ImportLinks.of(targets);
  return new Promise((resolve) => {
    const owner = Scope.root('import wizard');
    owner.build(() => new ImportFlow(table, targets, links, options, inserts, auto ?? null, auto === undefined,
      (report) => {
        owner.dispose();
        resolve(report);
      }).open());
  });
}

async function autoNumbered(table: DomainTable): Promise<string | null | undefined> {
  const at = table.address.indexOf('.');
  const manifest = await grok.dapi.domains.schema(table.address.slice(0, at)).manifest().catch(() => null);
  if (manifest === null)
    return undefined;
  const columns = manifest.tables[table.address.slice(at + 1)]?.columns ?? {};
  return table.info.businessKey.find((k) => columns[k]?.['autoNumber'] != null) ?? null;
}

class ImportFlow {
  private readonly _frame = signal<DG.DataFrame | null>(null);
  private readonly _rows = signal<Mapping[]>([]);
  private readonly _choices = new Map<string, ChoiceInput>();
  private readonly _recheck = signal(0);
  private readonly _report = div([], 'u2-domain-import-report');
  // absent where `support.batch` refuses the option
  private _mode?: ChoiceInput;
  private _allOrNothing?: BoolInput;
  private _errorOnDuplicate?: BoolInput;
  private _writeBack?: BoolInput;
  private _wizard!: Wizard;
  private _committed = false;
  private _result: Report | null = null;
  private _settled = false;
  // the preview posts again only when the key (the mapping, the mode and the flags) changed
  private _checked: {key: string, report: Promise<DomainBatchValidationLike | DomainBatchReportLike>} | undefined;
  // once per source and mapping: the importer can be slow, and the preview and the commit read the same plan
  private _planned: {frame: DG.DataFrame, key: string, plan: Promise<Planned[]>} | undefined;
  // an older dry run answering last must not paint over a newer one
  private _previewGen = 0;
  // '' while the dry run runs, a reason where all-or-nothing would abort
  private readonly _previewBlock = signal<string | null>(null);
  private readonly _importText = signal('Import');
  private _refusal: string[] = [];
  private readonly _aborted = signal(false);
  private _linked: {frame: DG.DataFrame, key: string, matches: Promise<Map<string, LinkMatch[]>>} | undefined;
  private readonly _linkStatus = signal<{text: string, ok: boolean}[]>([]);

  constructor(private readonly _table: DomainTable, private readonly _targets: IProperty[],
    private readonly _links: Map<string, ImportLink>, private readonly _options: DomainImportOptions,
    private readonly _inserts: boolean, private readonly _auto: string | null, private readonly _autoUnread: boolean,
    private readonly _resolve: (report: DomainBatchReportLike | null) => void) {}

  open(): void {
    this._wizard = new Wizard({
      steps: [
        // the source step gates on the source alone: the mapping it would gate on does not exist yet
        {id: 'source', title: 'Source', content: () => this._sourceStep(),
          canProceed: () => this._frame.value === null ? 'Choose an open table or a file' : null},
        {id: 'mapping', title: 'Mapping', content: () => this._mappingStep(), canProceed: () => this._gate()},
        {id: 'preview', title: 'Preview', content: () => this._previewStep(), nextText: this._importText,
          canProceed: () => this._gate() ?? this._previewBlock.value},
        {id: 'report', title: 'Report', content: () => this._report, done: true,
          onActivate: () => void this._commit(),
          actions: [{text: 'BACK', run: () => this._again(), visible: this._aborted}]},
      ],
      onFinish: () => this._settle(this._result),
      onCancel: () => this._settle(null),
    });
    this._wizard.root.classList.add('u2-domain-import');
    this._wizard.openInDialog(`Import ${nounOf(this._table.info)}`, {width: 720, height: 560});
  }

  private _gate(): string | null {
    return this._problems()[0] ?? null;
  }

  private _settle(report: DomainBatchReportLike | null): void {
    if (this._settled)
      return;
    this._settled = true;
    for (const choice of this._choices.values())
      choice.dispose();
    this._choices.clear();
    this._wizard.dispose();
    this._resolve(report);
  }

  private _sourceStep(): HTMLElement {
    const picker = fromDartInput<DG.DataFrame | null>(ui.input.table('Data'), 'source');
    picker.root.title = 'Choose an open table or a file';
    const open = grok.shell.tables;
    if (this._options.source !== undefined)
      picker.value = this._options.source;
    else if (open.length === 1)
      picker.value = open[0];
    const summary = span('', 'u2-domain-import-picked u2-form-note');
    const batch = this._table.table.support.batch;
    if (batch.upsert && this._table.info.businessKey.length > 0) {
      const upsert = {value: 'upsert', label: 'Add or update'};
      this._mode = new ChoiceInput({label: 'Mode', name: 'mode', value: this._inserts ? 'insert' : 'upsert',
        nullable: false, items: this._inserts ? [{value: 'insert', label: 'Add new rows'}, upsert] : [upsert]});
    }
    if (batch.partial) {
      this._allOrNothing = new BoolInput({label: 'All or nothing', value: true, name: 'allOrNothing',
        tooltipText: 'Any invalid row aborts the whole import'});
    }
    const duplicates = `A row whose ${this._keyHeader()} the table already has is an error instead of a skip`;
    if (batch.skipDuplicates) {
      this._errorOnDuplicate = new BoolInput({label: `Refuse ${nounOf(this._table.info)} that exist`,
        value: false, name: 'errorOnDuplicate', tooltipText: duplicates});
    }
    const column = this._keyColumn();
    if (column !== null)
      this._writeBack = new BoolInput({label: `Write ${this._keyHeader()} back`, value: true, name: 'writeBack'});
    const form = new Form().addAll([picker, this._mode, this._allOrNothing, this._errorOnDuplicate, this._writeBack]
      .filter((input) => input !== undefined));
    picker.root.parentElement!.insertBefore(summary, picker.root.nextSibling);
    form.effect(() => {
      const frame = picker.value;
      summary.textContent = frame === null ? '' :
        `${plural(frame.rowCount, 'row', 'rows')}, ${plural(frame.columns.names().length, 'column', 'columns')}`;
      if (this._writeBack !== undefined) {
        this._writeBack.visible = frame !== null;
        this._writeBack.tooltipText = frame === null ? '' :
          `Adds a "${column}" column to ${frame.name} with each row's ${this._keyHeader()}, as imported or as found`;
      }
      this._frame.value = frame;
      this._remap(frame);
      this._touch();
    });
    form.effect(() => {
      const insert = this._mode?.value !== 'upsert';
      if (this._errorOnDuplicate !== undefined) {
        this._errorOnDuplicate.disabled = insert ? null : 'Only when adding new rows';
        this._errorOnDuplicate.tooltipText = duplicates;
      }
      this._touch();
    });
    return divV([form.root], 'u2-domain-import-source');
  }

  private _mappingStep(): HTMLElement {
    const table = new DataTable<Mapping>({
      columns: [
        {name: 'column', header: 'Source column', width: 'minmax(120px, 1fr)'},
        {name: 'sample', header: 'Sample', render: (item) => this._sample(item.column)},
        {name: 'target', header: 'Import as', width: 'minmax(200px, 1fr)',
          render: (item) => this._choices.get(item.column)!.root},
      ],
      keyOf: (item) => item.column,
      rowHeight: 32,
      items: this._rows,
    });
    table.root.classList.add('u2-domain-import-mapping');
    const links = divV([]);
    links.dataset.u2Part = 'import-links';
    const hints = divV([]);
    const step = divV([table.root, links, hints, statusLine(signal(null), computed(() => this._mappingStatus()))],
      'u2-domain-import-mapping-step');
    const control = new Control(step);
    control.effect(() => {
      this._recheck.value;
      const frame = this._frame.value;
      if (frame !== null && this._links.size > 0)
        void this._checkLinks(frame);
    });
    control.effect(() => {
      this._recheck.value;
      const mapped = this._mapping();
      hints.replaceChildren(...this._targets.filter((p) => mapped.has(p.name!))
        .flatMap((p) => ImportLinks.idHint(p) ?? []).map((line) => divH([dot('accent'), span(line)])));
    });
    control.effect(() => links.replaceChildren(...this._linkStatus.value.map((line) =>
      divH([dot(line.ok ? 'success' : 'error'), span(line.text)]))));
    return step;
  }

  private async _checkLinks(frame: DG.DataFrame): Promise<void> {
    const mapping = JSON.stringify([...this._mapping()]);
    const matches = await this._resolveLinks(frame, [...this._mapping()]).catch((e) => e as Error);
    if (this._frame.peek() !== frame || mapping !== JSON.stringify([...this._mapping()]))
      return;
    const byLink = this._linkMapping([...this._mapping()]);
    const status = await Promise.all([...byLink.keys()].map(async (column) => {
      const link = this._links.get(column)!;
      const by = byLink.get(column)!.map(([field]) => DomainTable.captionOf(field.prop)).join(' + ');
      const title = `${DomainTable.captionOf(link.prop)} by ${by}`;
      if (matches instanceof Error)
        return {text: `${title}: ${matches.message}`, ok: false};
      const named = matches.get(column)!.filter((m) => m.id !== null || m.error !== null);
      const hit = named.filter((m) => m.id !== null).length;
      const found = `${hit === named.length ? 'all ' : `${hit} of `}${plural(named.length, 'row', 'rows')} found`;
      const byKey = hit === 0 && named.length > 0 ? await this._foundByKey(frame, link, byLink.get(column)!) : null;
      return {text: `${title}: ${found}${byKey ?? ''}`, ok: hit === named.length && hit > 0};
    }));
    if (this._frame.peek() === frame && mapping === JSON.stringify([...this._mapping()]))
      this._linkStatus.value = status;
  }

  private async _foundByKey(frame: DG.DataFrame, link: ImportLink, fields: [LinkField, string][]):
    Promise<string | null> {
    const key = link.target.info.businessKey;
    const field = key.length === 1 && fields.length === 1 && fields[0][0].path !== key[0] ?
      link.fields.find((f) => f.path === key[0]) : undefined;
    if (field === undefined)
      return null;
    const column = frame.columns.byName(fields[0][1]);
    const values = Array.from({length: frame.rowCount}, (_, i) => [ImportFlow._valueOrNull(column, i)]);
    const hit = (await ImportLinks.resolve(link, [field], values).catch(() => [])).filter((m) => m.id !== null).length;
    return hit === 0 ? null : `; as ${DomainTable.captionOf(link.prop)} › ${DomainTable.captionOf(field.prop)}, ` +
      `${plural(hit, 'row', 'rows')} would be found`;
  }

  // per link, the fields mapped to it and the file column each comes from
  private _linkMapping(mapping: [string, string][]): Map<string, [LinkField, string][]> {
    const out = new Map<string, [LinkField, string][]>();
    for (const [target, column] of mapping) {
      const [name, field] = target.split('.');
      const prop = field === undefined ? undefined : this._links.get(name)?.fields.find((f) => f.prop.name === field);
      if (prop !== undefined)
        out.set(name, [...out.get(name) ?? [], [prop, column]]);
    }
    return out;
  }

  private _resolveLinks(frame: DG.DataFrame, mapping: [string, string][]): Promise<Map<string, LinkMatch[]>> {
    const byLink = this._linkMapping(mapping);
    const key = JSON.stringify([...byLink].map(([name, fields]) => [name, fields.map(([f, c]) => [f.path, c])]));
    if (this._linked?.frame !== frame || this._linked.key !== key) {
      const matches = Promise.all([...byLink].map(async ([name, fields]): Promise<[string, LinkMatch[]]> => {
        const columns = fields.map(([, column]) => frame.columns.byName(column));
        const values = Array.from({length: frame.rowCount},
          (_, i) => columns.map((c) => ImportFlow._valueOrNull(c, i)));
        return [name, await ImportLinks.resolve(this._links.get(name)!, fields.map(([f]) => f), values)];
      })).then((entries) => new Map(entries));
      this._linked = {frame, key, matches};
      matches.catch(() => {
        if (this._linked?.matches === matches)
          this._linked = undefined;
      });
    }
    return this._linked.matches;
  }

  private _sample(name: string): HTMLElement {
    const column = this._frame.peek()?.columns.byName(name);
    const values: string[] = [];
    for (let i = 0; column && i < column.length && values.length < 3; i++) {
      if (!column.isNone(i))
        values.push(text(column.get(i)));
    }
    return span(values.join(', '), 'u2-domain-import-text');
  }

  private _mappingStatus(): string {
    this._recheck.value;
    const mapped = this._mapping();
    const imported = [...this._choices.values()].filter((c) => this._claimants(c.peek()).length > 0).length;
    const required = this._required()
      .map((name) => `${ImportFlow._caption(this, name)} ${ImportFlow._covers(mapped, name) ? '✓' : '✗'}`);
    return [...(required.length > 0 ? [`Required: ${required.join(', ')}`] : []),
      ...[...mapped.keys()].filter((t) => this._claimants(t).length > 1)
        .map((name) => `${ImportFlow._caption(this, name)} is mapped twice`),
      `${imported} of ${plural(this._choices.size, 'column', 'columns')} imported`].join(' · ');
  }

  private _previewStep(): HTMLElement {
    const host = div([], 'u2-domain-import-preview');
    const control = new Control(host);
    control.effect(() => {
      const frame = this._frame.value;
      const problems = this._problems();
      // the step is built once and kept: while the user is on another one, a keystroke there is
      // not a reason to post
      if (this._wizard.currentStep.value !== 'preview')
        return;
      host.replaceChildren(...problems.map((p) => span(p, 'u2-domain-import-problem')));
      if (frame === null || problems.length > 0)
        return;
      void this._preview(control, host, frame);
    });
    return host;
  }

  private async _preview(owner: Control, host: HTMLElement, frame: DG.DataFrame): Promise<void> {
    const mapping = [...this._mapping()];
    const options = this._batchOptions();
    const gen = ++this._previewGen;
    this._previewBlock.value = '';
    this._importText.value = 'Import';
    const status = span(`Checking ${plural(frame.rowCount, 'row', 'rows')}…`, 'u2-domain-import-summary');
    host.replaceChildren(status);
    let plan: Planned[];
    let report: DomainBatchValidationLike | null = null;
    let refusal: string | null = null;
    let notes: string[] = [];
    try {
      plan = await this._plan(frame, mapping,
        (done) => status.textContent = `Checking rows: ${done} of ${frame.rowCount}`);
      const sent = ImportFlow._sent(plan);
      if (this._table.table.validate !== undefined) {
        const key = JSON.stringify([mapping, options]);
        const scanned = sent.slice(0, PREVIEW_ROWS);
        const payloads = scanned.map((i) => plan[i].payload!);
        if (this._checked?.key !== key)
          this._checked = {key, report: this._table.table.validate(payloads, options)};
        const answer = await this._checked.report;
        if (ImportFlow._refused(answer)) {
          refusal = ImportFlow._refusalText(answer);
          const judged = await this._attributed(answer.rows, payloads);
          notes = ImportFlow._refusalNotes(judged.notes, scanned);
          report = ImportFlow._predicted(judged.rows, scanned.length);
        } else
          report = answer;
      }
      if (gen !== this._previewGen)
        return;
    } catch (e) {
      if (gen !== this._previewGen)
        return;
      // the commit is the authority: a preview that cannot reach the server says so and lets NEXT
      // through
      this._checked = undefined;
      this._previewBlock.value = null;
      host.replaceChildren(span(DomainErrors.message(e), 'u2-domain-import-problem'));
      return;
    }
    if (owner.scope.isDisposed)
      return;
    const sent = ImportFlow._sent(plan);
    const kept = plan.filter((p) => p.payload === null);
    const keptErrors = kept.filter((p) => p.verdict.error !== undefined).length;
    const errors = (report?.errorCount ?? 0) + keptErrors;
    const keys = this._writeBack?.peek() !== true ? 0 :
      kept.filter((p) => p.verdict.existing !== undefined).length +
      (report?.rows.filter((r) => r.predicted === 'skip' && r.existingId !== undefined).length ?? 0);
    const valid = report === null ? sent.length : sent.length > PREVIEW_ROWS ? null :
      report.willInsert + report.willUpdate;
    const noun = (n: number) => plural(n, nounOf(this._table.info, true), nounOf(this._table.info));
    const aborts = options.allOrNothing !== false && errors > 0;
    this._importText.value = refusal !== null ? 'Import' :
      aborts ? `Nothing to import while ${plural(errors, 'row has', 'rows have')} errors` :
        valid === 0 && keys > 0 ? `Write back ${plural(keys, 'key', 'keys')}` : `Import ${noun(valid ?? sent.length)}`;
    this._previewBlock.value = refusal !== null ? (notes.length > 0 ? 'Nothing would be imported until that is fixed.' :
      `The server refuses the whole batch — ${refusal}. Nothing would be imported until that is fixed.`) :
      plan.length === 0 ? 'Nothing would be imported — the file has no rows.' :
        aborts ? this._abortReason(errors, valid) : valid === 0 && keys === 0 ? ImportFlow._nothing(errors) : null;
    const judged = report !== null || this._table.importer?.prepare !== undefined;
    const rows = this._previewRows(plan, sent, report).slice(0, PREVIEW_SHOWN);
    host.replaceChildren(this._rowsTable(owner, 'preview', rows, judged ? this._previewColumns(mapping) :
      mapping.map(([target]) => this._valueColumn(target))));
    const showing = rows.length < plan.length ?
      `Showing the first ${rows.length.toLocaleString()} of ${plan.length.toLocaleString()} rows` : '';
    if (report === null) {
      host.append(span(`${showing === '' ? 'Rows' : `${showing} — rows`} are checked when imported.`,
        'u2-domain-import-summary'));
      return;
    }
    const checked = sent.length > PREVIEW_ROWS ? `Checked the first ${PREVIEW_ROWS.toLocaleString()} of ` +
      `${sent.length.toLocaleString()} rows — the rest are checked on import:` : '';
    host.prepend(ImportFlow._headline(ImportFlow._counts([[report.willInsert, 'to add', 'success'],
      [report.willUpdate, 'to update', 'success'], [report.willSkip + kept.length - keptErrors, 'to skip', 'default'],
      [errors, errors === 1 ? 'error' : 'errors', 'error']]), refusal === null ? 'Nothing to import' :
      'Refused by the server'));
    if (checked !== '')
      host.prepend(span(checked, 'u2-domain-import-summary'));
    if (showing !== '')
      host.append(span(`${showing}.`, 'u2-domain-import-summary'));
    host.append(...notes.map((note) => span(note, 'u2-domain-import-problem')));
  }

  private static _nothing(errors: number): string {
    return `Nothing would be imported — ${errors > 0 ? `${plural(errors, 'row has errors', 'rows have errors')}, ` +
      'and no row is valid' : 'no row would be added or updated'}.`;
  }

  private static _refused(answer: DomainBatchValidationLike | DomainBatchReportLike): answer is DomainBatchReportLike {
    return !('validateOnly' in answer);
  }

  private static _refusalText(refusal: DomainBatchReportLike): string {
    const messages = [...new Set(refusal.rows.flatMap((r) => (r.errors ?? []).map((e) => e.message)))];
    return messages.length > 0 ? messages.join('; ') : refusal.error ?? 'no reason given';
  }

  private static _predicted(rows: Placed[], rowCount: number): DomainBatchValidationLike {
    const refused = rows.filter((r) => r.errors?.length);
    return {validateOnly: true, rowCount, willInsert: 0, willUpdate: 0, willSkip: 0, errorCount: refused.length,
      rows: refused.map((r) => ({index: r.index, predicted: 'error', errors: r.errors}))};
  }

  private async _attributed(rows: DomainBatchRowLike[], payloads: Record<string, unknown>[]):
    Promise<{rows: Placed[], notes: {message: string, at: number[]}[]}> {
    const out = new Map<number, Placed>();
    const notes: {message: string, at: number[]}[] = [];
    const at = (index: number, from?: DomainBatchRowLike): Placed => {
      const row = out.get(index) ?? {...from, index, id: from?.id ?? null, errors: []};
      out.set(index, row);
      return row;
    };
    for (const row of rows) {
      if (row.index !== undefined) {
        at(row.index, row).errors!.push(...row.errors ?? []);
        continue;
      }
      for (const error of row.errors ?? []) {
        const breaking = error.code === 'check' ? await this._breaking(error.constraint, payloads) ?? [] : [];
        notes.push({message: this._describe([error]), at: breaking});
        for (const index of breaking)
          at(index).errors!.push(error);
      }
    }
    return {rows: [...out.values()].sort((a, b) => a.index - b.index), notes};
  }

  private async _breaking(name: string | undefined, payloads: Record<string, unknown>[]): Promise<number[] | null> {
    const admitted = await this._table.constraints.find((c) => c.name === name)?.admits(payloads);
    return admitted == null ? null : payloads.flatMap((_, i) => admitted[i] ? [] : [i]);
  }

  private static _refusalNotes(notes: {message: string, at: number[]}[], sent: number[]): string[] {
    return notes.map(({message, at}) => {
      const shown = at.slice(0, 10).map((i) => (sent[i] + 1).toLocaleString()).join(', ');
      const more = at.length > 10 ? ` and ${(at.length - 10).toLocaleString()} more` : '';
      return `The server refuses the whole batch — ${message}` + (at.length === 0 ? '; it does not say which row.' :
        ` (${at.length === 1 ? 'row' : 'rows'} ${shown}${more}).`);
    });
  }

  private static _headline(items: [string, BadgeVariant][], none: string): HTMLElement {
    return div(items.length === 0 ? [span(none, 'u2-domain-import-count')] : items.map(([words, variant]) =>
      div([dot(variant), span(words)], 'u2-domain-import-count')), 'u2-domain-import-outcome');
  }

  private static _counts(counts: [number, string, BadgeVariant][]): [string, BadgeVariant][] {
    return counts.filter(([n]) => n > 0).map(([n, words, variant]) => [`${n.toLocaleString()} ${words}`, variant]);
  }

  // `valid` is null where the dry run did not see every row
  private _abortReason(errors: number, valid: number | null): string {
    const lead = plural(errors, 'row has errors', 'rows have errors');
    if (this._allOrNothing === undefined)
      return `${lead}, and this import is all or nothing.`;
    if (valid === 0)
      return `${lead}, and no row is valid.`;
    return `${lead}. Uncheck "All or nothing" on Source to import the ` +
      `${valid === null ? 'valid rows' : plural(valid, 'valid row', 'valid rows')} and skip the rest.`;
  }

  private _previewRows(plan: Planned[], sent: number[], report: DomainBatchValidationLike | null): DomainImportRow[] {
    const rows = plan.map((p, i) => p.payload === null ? this._kept(i, p) :
      this._row(i, p, report === null ? 'insert' : '', p.verdict.note ?? '', []));
    for (const dry of report?.rows ?? []) {
      const at = sent[dry.index];
      const errors = dry.errors ?? [];
      rows[at] = this._row(at, plan[at], dry.predicted, errors.length > 0 ? this._describe(errors) :
        plan[at].verdict.note ?? (dry.predicted === 'skip' ? this._duplicate() : ''), errors);
    }
    return rows;
  }

  private _kept(index: number, planned: Planned): DomainImportRow {
    const verdict = planned.verdict;
    return verdict.error !== undefined ? this._row(index, planned, 'error', verdict.error, []) :
      {...this._row(index, planned, 'skip', verdict.note ?? '', []), id: verdict.existing ?? null};
  }

  private _row(index: number, planned: Planned, kind: string, note: string, errors: DomainImportError[]):
    DomainImportRow {
    const typed = Object.fromEntries(this._targets.filter((p) => p.semType === 'User' && p.name! in planned.cells)
      .map((p) => [p.name!, planned.cells[p.name!]]));
    return {index, values: {...planned.cells, ...planned.payload, ...typed}, kind, label: planned.verdict.label ?? '',
      note, errors, id: null, key: ''};
  }

  private _rowsTable(owner: Control, step: DomainImportStep, rows: DomainImportRow[],
    defaults: DomainImportColumn[]): HTMLElement {
    const keyless = rows.every((row) => row.key === '');
    const columns = (this._table.importer?.columns?.(defaults, step) ?? defaults)
      .filter((c) => c.name !== '~key' || !keyless);
    const byIndex = new Map(rows.map((row) => [row.index, row]));
    const table = owner.build(() => new DataTable<DomainImportRow>({
      rowHeight: this._table.importer?.rowHeight ?? 24,
      keyOf: (row) => String(row.index),
      fitColumns: true,
      columns: columns.map((c) => ({name: c.name, header: c.header, width: c.width, render: (row) => {
        const content = c.render(row);
        return typeof content === 'string' ? span(content, 'u2-domain-import-text') : content;
      }})),
      cellState: {
        isChanged: () => false,
        errorOf: (key, column) => {
          const message = byIndex.get(Number(key))?.errors.find((e) => e.column === column)?.message;
          return message === undefined ? null : {message, kind: 'error'};
        },
      },
    }));
    table.setItems(rows);
    table.root.classList.add('u2-domain-import-rows', `u2-domain-import-${step}-rows`);
    return table.root;
  }

  private _previewColumns(mapping: [string, string][]): DomainImportColumn[] {
    return [ImportFlow._numberColumn(), ImportFlow._resultColumn(VERDICTS), ImportFlow._noteColumn(),
      ...mapping.map(([target]) => this._valueColumn(target))];
  }

  private _reportColumns(saved: Map<string, RowView>): DomainImportColumn[] {
    return [ImportFlow._numberColumn(), ImportFlow._resultColumn(OUTCOMES),
      {name: '~key', header: this._keyHeader(), width: 'minmax(120px, auto)',
        render: (row) => this._keyLink(row, saved)},
      ImportFlow._noteColumn()];
  }

  private _valueColumn(target: string): DomainImportColumn {
    const prop = this._table.properties.find((p) => p.name === target);
    return {name: target, header: ImportFlow._caption(this, target),
      width: DomainDataTable.width(prop, target === this._table.info.nameColumn),
      render: (row) => text(row.values[target])};
  }

  private static _numberColumn(): DomainImportColumn {
    return {name: '~row', header: '#', width: '44px', render: (row) => String(row.index + 1)};
  }

  private static _noteColumn(): DomainImportColumn {
    return {name: '~note', header: 'Note', width: 'minmax(160px, 2fr)', render: (row) => row.note};
  }

  private static _resultColumn(kinds: Record<string, {label: string, variant: BadgeVariant}>): DomainImportColumn {
    return {name: '~result', header: 'Result', width: '88px', render: (row) => {
      const kind = kinds[row.kind];
      const label = kind?.label ?? row.label;
      if (!label)
        return '';
      const pill = badge(label, {variant: kind?.variant ?? 'default'});
      return kind?.variant === 'error' ? div([dot('error'), pill], 'u2-domain-import-result') : pill;
    }};
  }

  private _keyLink(row: DomainImportRow, saved: Map<string, RowView>): HTMLElement | string {
    const target = row.id === null ? undefined : saved.get(row.id);
    const name = this._table.info.nameColumn;
    return target === undefined ? (name == null ? '' : text(row.values[name])) : link(row.key, () => {
      this._settle(this._result);
      this._table.open(target);
    });
  }

  private _describe(errors: DomainImportError[]): string {
    return errors.map((e) => {
      const caption = e.column ? ImportFlow._caption(this, e.column) : '';
      return this._table.importer?.describe?.(e, caption) ?? (e.code === 'unique' && caption ?
        `${caption} is already used` : caption ? `${caption}: ${e.message}` : e.message);
    }).join('; ');
  }

  // the table's validators too: an import is held to the rules a form holds a row to
  private _plan(frame: DG.DataFrame, mapping: [string, string][], progress: (done: number) => void):
    Promise<Planned[]> {
    const key = JSON.stringify(mapping);
    if (this._planned?.frame !== frame || this._planned.key !== key) {
      const plan = this._planRows(frame, mapping, progress);
      this._planned = {frame, key, plan};
      plan.catch(() => {
        if (this._planned?.plan === plan)
          this._planned = undefined;
      });
    }
    return this._planned.plan;
  }

  private async _planRows(frame: DG.DataFrame, mapping: [string, string][], progress: (done: number) => void):
    Promise<Planned[]> {
    const cells = this._payload(frame, mapping, frame.rowCount);
    // what is sent: a link's id in place of the fields that found it; a row naming no record leaves the link alone
    const rows = cells.map((row) => Object.fromEntries(Object.entries(row).filter(([key]) => !key.includes('.'))));
    const missed: string[][] = rows.map(() => []);
    for (const [column, matches] of await this._resolveLinks(frame, mapping)) {
      for (const [i, row] of rows.entries()) {
        if (matches[i].id !== null)
          row[column] = matches[i].id;
        if (matches[i].error !== null)
          missed[i].push(matches[i].error!);
      }
    }
    await this._namedUsers(rows);
    const prepare = this._table.importer?.prepare;
    const verdicts = prepare === undefined ? [] : await prepare(rows, progress);
    return cells.map((cell, i): Planned => {
      const verdict = missed[i].length > 0 ?
        {error: [...missed[i], ...this._validate(rows[i]), ...this._valueProblems(rows[i])].join('; ')} :
        verdicts[i] ?? {};
      if (verdict.error !== undefined || verdict.existing !== undefined || verdict.sameAs !== undefined)
        return {cells: cell, payload: null, verdict};
      const payload = {...rows[i], ...verdict.values};
      const problems = this._validate(payload);
      return problems.length === 0 ? {cells: cell, payload, verdict} :
        {cells: cell, payload: null, verdict: {error: problems.join('; ')}};
    });
  }

  private async _namedUsers(rows: Record<string, unknown>[]): Promise<void> {
    const columns = this._targets.filter((p) => p.semType === 'User');
    const named = [...new Set(rows.flatMap((row) => columns.flatMap((p) => {
      const said = ImportLinks.idProblem(p, row[p.name!]) === null ? '' : String(row[p.name!]).trim();
      return said === '' ? [] : [said];
    })))];
    const users = new Map<string, DG.User>();
    const filters = ImportFlow._userFilters(named);
    for (let at = 0; at < filters.length; at += USER_READS) {
      await Promise.all(filters.slice(at, at + USER_READS).map(async (filter) => {
        for (const user of await grok.dapi.users.filter(filter).list())
          users.set(user.id, user);
      }));
    }
    const ids = new Map<string, string | null>();
    for (const {id, friendlyName, firstName, lastName} of users.values()) {
      const names = [friendlyName, firstName, lastName, `${firstName ?? ''} ${lastName ?? ''}`.trim()];
      for (const key of new Set(names.filter((n) => !!n).map((n) => n.toLowerCase())))
        ids.set(key, ids.has(key) ? null : id);
    }
    for (const user of [...users.values()].filter((u) => !!u.login))
      ids.set(user.login.toLowerCase(), user.id);
    for (const row of rows) {
      for (const p of columns) {
        const id = ImportLinks.idProblem(p, row[p.name!]) === null ? null :
          ids.get(String(row[p.name!]).trim().toLowerCase());
        if (id)
          row[p.name!] = id;
      }
    }
  }

  private static _userClauses(v: string): string[] {
    const q = (s: string) => JSON.stringify(s);
    const words = v.split(/\s+/);
    return [...new Set([v, v.toLowerCase()])].map((login) => `login = ${q(login)}`).concat(
      `firstName = ${q(v)}`, `lastName = ${q(v)}`, ...words.slice(1).map((_, i) =>
        `(firstName = ${q(words.slice(0, i + 1).join(' '))} and lastName = ${q(words.slice(i + 1).join(' '))})`));
  }

  private static _userFilters(named: string[]): string[] {
    const size = (s: string) => encodeURIComponent(encodeURIComponent(s)).length;
    const or = size(' or ');
    const filters: string[][] = [];
    let bytes = Infinity;
    for (const v of named) {
      const clauses = ImportFlow._userClauses(v);
      const more = clauses.reduce((n, c) => n + size(c) + or, 0);
      if (more > USER_FILTER_BYTES)
        continue;
      if (bytes + more > USER_FILTER_BYTES) {
        filters.push([]);
        bytes = 0;
      }
      filters[filters.length - 1].push(...clauses);
      bytes += more;
    }
    return filters.map((clauses) => clauses.join(' or '));
  }

  private _validate(row: Record<string, unknown>): string[] {
    const validators = this._table.validators;
    const problems = this._targets.map((prop) => [prop.name!, ImportLinks.idProblem(prop, row[prop.name!])]);
    if (this._auto !== null && this._auto in row && isEmpty(row[this._auto]))
      problems.push([this._auto, 'A new row needs Add new rows mode']);
    for (const column of validators.columns)
      problems.push([column, validators.check(column, row[column], row as RowView)]);
    return problems.filter(([, problem]) => problem !== null)
      .map(([column, problem]) => `${ImportFlow._caption(this, column!)}: ${problem}`);
  }

  private _valueProblems(row: Record<string, unknown>): string[] {
    return this._targets.filter((prop) => prop.name! in row && !ImportLinks.isLink(prop)).flatMap((prop) => {
      const value = row[prop.name!];
      const numeric = ['int', 'bigint', 'double', 'float'].includes(prop.propertyType ?? prop.type ?? '');
      const problem = numeric && typeof value === 'string' && !Number.isFinite(Number(value)) ?
        `"${value}" is not a number` : MemoryEditState.problemOf(prop, value);
      return problem === null ? [] : [`${ImportFlow._caption(this, prop.name!)}: ${problem}`];
    });
  }

  // a sent row's place in the batch is its index here
  private static _sent(plan: Planned[]): number[] {
    return plan.flatMap((p, i) => p.payload === null ? [] : [i]);
  }

  // what is posted, and the only thing that is: a skipped or renamed source column never reaches the server
  private _payload(frame: DG.DataFrame, mapping: [string, string][], count: number):
    Record<string, unknown>[] {
    const rows: Record<string, unknown>[] = Array.from({length: count}, () => ({}));
    for (const [target, name] of mapping) {
      const column = frame.columns.byName(name);
      const key = this._table.info.businessKey.includes(target) ? this._targets.find((p) => p.name === target &&
        p.format && ['int', 'bigint'].includes(p.propertyType ?? p.type ?? '')) : undefined;
      for (let i = 0; i < count; i++) {
        const value = ImportFlow._valueOrNull(column, i);
        rows[i][target] = key !== undefined && typeof value === 'string' ? DomainTable.parseKey(key, value) ?? value :
          value;
      }
    }
    return rows;
  }

  private static _valueOrNull(column: DG.Column, i: number): unknown {
    return column.isNone(i) ? null : column.get(i);
  }

  // an option the form does not offer is left to the server's default
  private _batchOptions(): DomainBatchOptionsLike {
    const mode = this._mode?.peek() === 'upsert' ? 'upsert' : 'insert';
    const options: DomainBatchOptionsLike = {mode};
    if (this._allOrNothing !== undefined)
      options.allOrNothing = this._allOrNothing.peek();
    if (this._errorOnDuplicate !== undefined)
      options.errorOnDuplicate = mode === 'insert' && this._errorOnDuplicate.peek();
    return options;
  }

  private _touch(): void {
    this._recheck.value = this._recheck.peek() + 1;
  }

  // the choices live here, not in the mapping table's cells: the table recycles its rows, and a choice in a
  // cell would be rebuilt, losing its value, on every scroll
  private _remap(frame: DG.DataFrame | null): void {
    this._checked = undefined;
    this._planned = undefined;
    this._linked = undefined;
    for (const choice of this._choices.values())
      choice.dispose();
    this._choices.clear();
    if (frame === null) {
      this._rows.value = [];
      return;
    }
    // the asterisk the forms use, in the option itself: a picker's own list is all the user sees
    const linked = [...this._links.values()].flatMap((link) => link.fields.map((f): ChoiceItem => {
      const title = DomainTable.captionOf(link.prop);
      return {value: `${link.prop.name}.${f.prop.name}`, label: `${title} › ${DomainTable.captionOf(f.prop)}`,
        group: `${title}: find by`};
    }));
    const items: ChoiceItem[] = [{value: SKIP, label: SKIP_LABEL}, ...this._targets.map((p) =>
      ({value: p.name!, label: `${DomainTable.captionOf(p)}${p.nullable === false ? ' *' : ''}`})), ...linked];
    const columns = frame.columns.names();
    for (const column of columns) {
      // owned by hand, not by a scope: the first remap runs inside the wizard's own constructor
      // (a source handed in as an option), where there is no wizard to own anything yet
      const choice = new ChoiceInput({items, value: this._match(frame.columns.byName(column)), inline: true,
        name: column, nullable: false});
      // required as the forms say it: a target the table will not take empty
      choice.effect(() => {
        const target = choice.value;
        const name = (target ?? '').split('.')[0];
        choice.required = this._targets.some((p) => p.name === name && p.nullable === false);
        choice.root.classList.toggle('u2-domain-import-skipped', target === SKIP);
      });
      choice.addValidator((target) => {
        const other = this._claimants(target).find((c) => c !== column);
        return other === undefined ? null : `Also mapped from "${other}"`;
      });
      choice.signal.subscribe(() => {
        this._touch();
        for (const each of this._choices.values())
          each.revalidate();
      });
      this._choices.set(column, choice);
    }
    this._rows.value = columns.map((column) => ({column}));
  }

  private _claimants(target: string | null): string[] {
    if (target === null || target === SKIP)
      return [];
    return [...this._choices].filter(([, choice]) => choice.peek() === target).map(([column]) => column);
  }

  private static _caption(flow: ImportFlow, column: string): string {
    const [name, field] = column.split('.');
    const prop = flow._table.properties.find((p) => p.name === name);
    const caption = prop === undefined ? name : DomainTable.captionOf(prop);
    const linked = flow._links.get(name)?.fields.find((f) => f.prop.name === field);
    return linked === undefined ? caption : `${caption} › ${DomainTable.captionOf(linked.prop)}`;
  }

  private _keyColumn(): string | null {
    const chosen = this._table.importer?.keyColumn;
    const key = this._table.info.businessKey;
    return chosen === false ? null : chosen ?? (key.length === 1 ? ImportFlow._caption(this, key[0]) : null);
  }

  private _duplicate(): string {
    return `${this._keyHeader()} already exists`;
  }

  private _keyHeader(): string {
    const key = this._table.info.businessKey;
    return key.length === 1 ? ImportFlow._caption(this, key[0]) :
      capital(nounOf(this._table.info, true));
  }

  // under the key header: the key's value as its format spells it (`DG-000060`), not the caption a name column gives
  private _keyOf(row: RowView): string {
    const key = this._table.info.businessKey;
    if (key.length !== 1)
      return this._table.renderer.caption(row);
    const value = row[key[0]];
    return typeof value === 'number' ?
      DomainDataTable.number(value, this._table.properties.find((p) => p.name === key[0])?.format) : text(value);
  }

  // "CAS number" is the cas column: case and non-alphanumerics drop, and an exact hit beats an earlier partial one
  private _match(column: DG.Column): string {
    const chosen = this._table.importer?.match?.(column, this._targets);
    if (chosen !== undefined && this._targets.some((p) => p.name === chosen))
      return chosen;
    const plain = (s: string) => s.toLowerCase().replace(/[^a-z0-9]/g, '');
    const want = plain(column.name);
    if (want === '')
      return SKIP;
    const names = (p: IProperty) => [plain(p.name ?? ''), plain(p.friendlyName ?? '')].filter((n) => n !== '');
    const hit = (test: (a: string) => boolean) => this._targets.find((p) => names(p).some(test))?.name;
    const typed = column.semType ? this._targets.filter((p) => p.semType === column.semType) : [];
    const found = hit((a) => a === want) ?? hit((a) => a.startsWith(want) || want.startsWith(a)) ??
      hit((a) => a.includes(want) || want.includes(a)) ?? (typed.length === 1 ? typed[0].name! : SKIP);
    const target = found === this._auto && this._mode?.peek() !== 'upsert' ? SKIP : found;
    const link = this._links.get(target);
    if (link === undefined || plain(link.prop.name!) === want || ImportFlow._ids(column))
      return target;
    const ref = names(link.prop).flatMap((n) => [n, n.replace(/id$/, '')]).filter((r) => r !== '' && want.includes(r))
      .sort((x, y) => y.length - x.length)[0];
    if (ref === undefined)
      return target;
    const rest = want.slice(want.indexOf(ref) + ref.length);
    const info = link.target.info;
    const field = (name?: string | null) => link.fields.find((f) => f.prop.name === name);
    const bare = want === ref && !/[^\x20-\x7e]/.test(column.name);
    const by = rest !== '' ? link.fields.find((f) => names(f.prop).includes(rest)) :
      bare ? field(info.nameColumn) ?? field(info.businessKey[0]) : undefined;
    return by === undefined ? target : `${target}.${by.prop.name}`;
  }

  private static _ids(column: DG.Column): boolean {
    for (let i = 0; i < column.length; i++) {
      if (!column.isNone(i) && !DomainAddress.ID.test(text(column.get(i))))
        return false;
    }
    return true;
  }

  private _mapping(): Map<string, string> {
    const out = new Map<string, string>();
    for (const [column, choice] of this._choices) {
      const target = choice.peek();
      if (target !== null && target !== SKIP && !out.has(target))
        out.set(target, column);
    }
    return out;
  }

  // in the order `domain_import_dialog.dart` reports them
  private _problems(): string[] {
    this._recheck.value;
    if (this._frame.value === null)
      return ['Choose an open table or a file'];
    const problems: string[] = [];
    const seen = this._mapping();
    for (const name of this._links.keys()) {
      if (seen.has(name) && [...seen.keys()].some((target) => target.startsWith(`${name}.`)))
        problems.push(`"${ImportFlow._caption(this, name)}" is mapped both itself and by its fields — keep one.`);
    }
    for (const target of seen.keys()) {
      const columns = this._claimants(target).map((c) => `"${c}"`);
      if (columns.length > 1) {
        problems.push(`Columns ${columns.slice(0, -1).join(', ')} and ${columns.at(-1)} are ` +
          `${columns.length === 2 ? 'both' : 'all'} mapped to "${ImportFlow._caption(this, target)}" — ` +
          `keep one, and set the rest to "${SKIP_LABEL}".`);
      }
    }
    if (seen.size === 0)
      problems.push('Map at least one column.');
    const upsert = this._mode?.value === 'upsert';
    if (!upsert && this._auto !== null && seen.has(this._auto)) {
      problems.push(`"${ImportFlow._caption(this, this._auto)}" is assigned on save — choose "Add or update", or set ` +
        `it to "${SKIP_LABEL}".`);
    }
    for (const name of this._required()) {
      if (ImportFlow._covers(seen, name))
        continue;
      const caption = ImportFlow._caption(this, name);
      const unread = upsert && this._autoUnread && !this._targets.some((p) => p.name === name);
      problems.push(unread ? `Could not read which key the server numbers: Add or update by ${caption} is off.` :
        upsert ? `Add or update matches rows by ${caption} — map a column to "${caption}".` :
          `Required column "${caption}" is not mapped — rows without a value will fail.`);
    }
    return problems;
  }

  private static _covers(mapping: Map<string, string>, name: string): boolean {
    return mapping.has(name) || [...mapping.keys()].some((target) => target.startsWith(`${name}.`));
  }

  private _required(): string[] {
    return this._mode?.value === 'upsert' ? this._table.info.businessKey : this._targets
      .filter((p) => p.nullable === false && (p.defaultValue ?? null) === null).map((p) => p.name!);
  }

  private async _commit(): Promise<void> {
    if (this._committed)
      return;
    this._committed = true;
    const frame = this._frame.peek()!;
    const mapping = [...this._mapping()];
    const options = this._batchOptions();
    this._report.replaceChildren(span(`Importing ${frame.rowCount} rows…`, 'u2-domain-import-summary'));
    let plan: Planned[];
    // reading the frame is part of the import: a column the picked frame no longer has belongs in
    // the report with everything else the server refuses
    try {
      plan = await this._plan(frame, mapping, () => {});
      const sent = ImportFlow._sent(plan);
      // a row the plan refuses is an error like the server's: all-or-nothing sends nothing then
      const refused = plan.some((p) => p.verdict.error !== undefined) && options.allOrNothing !== false;
      const payloads = sent.map((i) => plan[i].payload!);
      const answer: DomainBatchReportLike = refused ?
        {inserted: 0, updated: 0, skipped: 0, errorCount: 0, error: 'validation', rows: []} :
        await this._table.table.batch!(payloads, options);
      const judged = await this._attributed(answer.rows, payloads);
      this._refusal = ImportFlow._refusalNotes(judged.notes, sent);
      this._result = ImportFlow._merged({...answer, rows: judged.rows}, sent, plan);
    } catch (e) {
      this._report.replaceChildren(span(DomainErrors.message(e), 'u2-domain-import-problem'));
      this._wizard.setStepFailed('report', true);
      // a balloon only for a dialog closed before the end
      if (this._settled)
        notify.error(DomainErrors.message(e));
      return;
    }
    await this._render(this._result, plan, frame);
  }

  private static _merged(report: Report, sent: number[], plan: Planned[]): Report {
    const kept = plan.flatMap((p, i): Placed[] => p.payload !== null ? [] :
      [p.verdict.error !== undefined ? {index: i, id: null, status: 'error', errors: [{message: p.verdict.error}]} :
        {index: i, id: null, status: 'skipped', existingId: p.verdict.existing}]);
    const errors = kept.filter((r) => r.status === 'error').length;
    return {...report, skipped: (report.skipped ?? 0) + kept.length - errors,
      errorCount: report.errorCount === undefined ? undefined : report.errorCount + errors,
      rows: [...report.rows.map((r) => ({...r, index: sent[r.index]})), ...kept].sort((a, b) => a.index - b.index)};
  }

  private async _render(report: Report, plan: Planned[], frame: DG.DataFrame): Promise<void> {
    const failed = report.error !== undefined;
    // a warehouse refusal carries the failing rows and no totals
    const errorCount = report.errorCount ?? report.rows.filter((row) => row.errors?.length).length;
    const counts = ImportFlow._counts(ImportFlow._landed(report, errorCount));
    const aborted = errorCount === 0 && this._refusal.length > 0 ?
      'Import aborted — the server refused the whole batch; nothing was committed.' :
      `Import aborted — ${plural(errorCount, 'row has errors', 'rows have errors')}; nothing was committed.`;
    this._report.replaceChildren(ImportFlow._headline(failed ? [[aborted, 'error']] : counts, 'Nothing was imported'));
    this._report.classList.toggle('u2-domain-import-aborted', failed);
    this._wizard.setStepFailed('report', failed || errorCount > 0);
    if (failed)
      this._report.append(...this._refusal.map((note) => span(note, 'u2-domain-import-problem')));
    let itemized = report;
    let saved = new Map<string, RowView>();
    try {
      itemized = await this._itemized(report, plan);
    } catch (e) {
      this._report.append(span(`The rows the server did not itemize could not be read: ${DomainErrors.message(e)}`,
        'u2-domain-import-problem'));
    }
    const outcomes = this._outcomes(itemized, plan, failed);
    try {
      saved = await this._rowsBy('id', outcomes.map((o) => o.id).filter((id) => id !== null));
    } catch (e) {
      this._report.append(span(`The keys could not be read: ${DomainErrors.message(e)}`, 'u2-domain-import-problem'));
    }
    for (const row of outcomes)
      row.key = row.id === null || !saved.has(row.id) ? '' : this._keyOf(saved.get(row.id)!);
    const unitemized = failed ?
      Math.max(0, errorCount - report.rows.filter(ImportFlow._isError).length) :
      outcomes.filter((o) => o.kind === '').length;
    if (unitemized > 0) {
      this._report.append(span(`${plural(unitemized, 'row is', 'rows are')} not itemized by the server: ` +
        (failed ? 'their errors are not listed.' : 'their outcome and key are unknown.'), 'u2-domain-import-problem'));
    }
    // a dialog closed while the keys were read has no report to show; the write-back and the balloon still follow
    if (!this._wizard.scope.isDisposed)
      this._report.append(this._rowsTable(this._wizard, 'report', outcomes, this._reportColumns(saved)));
    const column = this._keyColumn();
    const keyed = outcomes.some((o) => o.key !== '');
    if (!failed && keyed && column !== null && this._writeBack?.peek()) {
      if (frame.columns.names().includes(column))
        frame.columns.remove(column);
      // column names ignore case: a source's own "sku" keeps its values beside the "Sku" written back
      const added = frame.columns.addNewString(frame.columns.getUnusedName(column)).init((i) => outcomes[i].key);
      this._report.append(span(`Added a "${added.name}" column to ${frame.name} with each row's ${this._keyHeader()}` +
        (unitemized > 0 ? `; it is empty for the ${plural(unitemized, 'row', 'rows')} not itemized.` : '.'),
        'u2-domain-import-summary'));
    }
    if (failed) {
      const refused = report.rows.filter((r) => ImportFlow._isError(r) && plan[r.index].payload !== null).length;
      const valid = ImportFlow._sent(plan).length - refused;
      if (this._allOrNothing?.peek() && valid > 0 && this._refusal.length === 0) {
        this._report.append(span(`Uncheck "All or nothing" to import the ${plural(valid, 'valid row', 'valid rows')} ` +
          'and skip the rest.', 'u2-domain-import-summary'));
      }
      this._aborted.value = true;
    }
    if (!this._settled)
      return;
    if (failed)
      notify.error(aborted);
    else
      notify.info(`Imported: ${counts.map(([words]) => words).join(', ')}`);
  }

  private static _isError(row: DomainBatchRowLike): boolean {
    return row.status === 'error' || (row.errors?.length ?? 0) > 0;
  }

  // a repeat takes its first row's id
  private _outcomes(report: Report, plan: Planned[], failed: boolean): DomainImportRow[] {
    const bySource = new Map(report.rows.map((r) => [r.index, r]));
    const outcomes = plan.map((p, i): DomainImportRow => {
      const row = bySource.get(i);
      const errors = row?.errors ?? [];
      if (row !== undefined && ImportFlow._isError(row))
        return this._row(i, p, 'error', p.payload === null ? p.verdict.error ?? '' : this._describe(errors), errors);
      if (row?.status === 'duplicate')
        return {...this._row(i, p, 'skip', this._duplicate(), errors), id: row.existingId ?? null};
      if (p.payload === null)
        return this._kept(i, p);
      if (failed || row === undefined) {
        return this._row(i, p, failed ? 'aborted' : '',
          failed ? 'Another row has errors' : 'Not itemized by the server', errors);
      }
      return {...this._row(i, p, row.status === 'inserted' ? 'insert' : 'update', p.verdict.note ?? '', errors),
        id: row.id};
    });
    plan.forEach((p, i) => {
      if (p.verdict.sameAs !== undefined)
        outcomes[i].id = outcomes[p.verdict.sameAs]?.id ?? null;
    });
    return outcomes;
  }

  private async _itemized(report: Report, plan: Planned[]): Promise<Report> {
    const listed = new Set(report.rows.map((r) => r.index));
    const missing = plan.flatMap((p, i) => p.payload === null || listed.has(i) ? [] : [i]);
    const key = this._table.info.businessKey;
    if (report.error !== undefined || missing.length === 0)
      return report;
    const sent = report.rows.filter((r) => plan[r.index].payload !== null);
    const keptErrors = plan.filter((p) => p.payload === null && p.verdict.error !== undefined).length;
    const keptSkips = plan.filter((p) => p.payload === null).length - keptErrors;
    const listedAs = (status: string) => sent.filter((r) => r.status === status).length;
    const left = ([['error', (report.errorCount ?? 0) - keptErrors], ['duplicate', (report.skipped ?? 0) - keptSkips],
      ['updated', report.updated ?? 0], ['merged', report.merged ?? 0], ['inserted', report.inserted ?? 0]] as const)
      .map(([status, total]) => [status, total - listedAs(status)] as const).filter(([, count]) => count > 0);
    if (left.length !== 1 || left[0][1] !== missing.length)
      return report;
    const status = left[0][0];
    if (status === 'error') {
      return {...report, rows: [...report.rows, ...missing.map((index) => ({index, id: null, status,
        errors: [{message: 'The server did not itemize this row\'s errors'}]}))]};
    }
    if (key.length !== 1)
      return report;
    const value = (index: number) => plan[index].payload![key[0]];
    const rows = await this._rowsBy(key[0], missing.map(value).filter((v) => v !== null && v !== undefined));
    const found = missing.flatMap((index) => {
      const row = rows.get(String(value(index)));
      const id = String(row?.id);
      return row === undefined ? [] : [{index, id, status, ...(status === 'duplicate' ? {existingId: id} : {})}];
    });
    return {...report, rows: [...report.rows, ...found]};
  }

  private async _rowsBy(column: string, values: unknown[]): Promise<Map<string, RowView>> {
    const rows = new Map<string, RowView>();
    const unique = [...new Set(values)];
    for (let at = 0; at < unique.length; at += KEY_CHUNK) {
      const chunk = unique.slice(at, at + KEY_CHUNK);
      const filter = {property: column, operator: '=', value: chunk};
      for (const row of await this._table.table.query({filter, limit: chunk.length}))
        rows.set(String(row[column]), row as RowView);
    }
    return rows;
  }

  // `merged` stands in for `updated` where the storage cannot tell the two apart
  private static _landed(report: DomainBatchReportLike, errors: number): [number, string, BadgeVariant][] {
    return [[report.inserted ?? 0, 'added', 'success'], [report.updated ?? 0, 'updated', 'success'],
      [report.merged ?? 0, 'merged', 'success'], [report.skipped ?? 0, 'skipped', 'default'],
      [errors, errors === 1 ? 'error' : 'errors', 'error']];
  }

  private _again(): void {
    this._committed = false;
    this._result = null;
    this._checked = undefined;
    this._aborted.value = false;
    this._wizard.setStepFailed('report', false);
    this._wizard.goTo('source');
  }
}

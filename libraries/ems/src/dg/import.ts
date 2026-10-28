import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import type * as DG from 'datagrok-api/dg';
import {signal, computed, Control, div, divV, span, plural, text, link, Wizard, DataTable, ChoiceInput,
  BoolInput, badge, dot, statusLine, Form, notify} from '@datagrok-libraries/u2';
import type {IProperty, ChoiceItem, BadgeVariant, RowView} from '@datagrok-libraries/u2';
import {capital, nounOf} from '../core/domain-backend.js';
import type {DomainBatchOptionsLike, DomainBatchReportLike, DomainBatchValidationLike} from '../core/domain-backend.js';
import {fromDartInput} from '@datagrok-libraries/u2/dg';
import {DomainTable} from './facade.js';
import {DomainDataTable} from './grid.js';
import {DomainErrors} from './errors.js';

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
   * `~result`, each mapped target by its name and `~note` in the preview; `~row`, `~result`, `~key` and
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

/** One row of the mapping step — the target is the choice input the flow keeps for that column. */
interface Mapping {column: string}

/** Opens the wizard over `table`; resolves to the server's report, or to null where it was
 * cancelled. Exported as `domains.import` — `import` cannot name a function. Only the mapped columns are posted,
 * under their target names, plus what the table's {@link DomainImporter} derives; the preview is the server's
 * dry run, the commit's own checks inside a rolled-back transaction. */
export function openImport(table: DomainTable, options: DomainImportOptions = {}):
  Promise<DomainBatchReportLike | null> {
  if (table.table.batch === undefined) {
    notify.error(`${table.address}: the backend does not support batch import`);
    return Promise.resolve(null);
  }
  if (!table.table.support.transaction) {
    notify.error(`${table.address}: the table does not accept writes`);
    return Promise.resolve(null);
  }
  // an import inserts: the draft policy, so a key typed once on insert is a target too
  const writable = table.access.columnPolicy().forDraft();
  const derived = table.importer?.derives ?? [];
  const targets = table.properties.filter((p) => writable.field(p.name!) === 'editable' && !derived.includes(p.name!));
  if (targets.length === 0) {
    notify.warning(`There is no column of ${table.address} you may write.`);
    return Promise.resolve(null);
  }
  return new Promise((resolve) => new ImportFlow(table, targets, options, resolve).open());
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
  private _result: DomainBatchReportLike | null = null;
  private _settled = false;
  // the preview posts again only when the key (the mapping, the mode and the flags) changed
  private _checked: {key: string, report: Promise<DomainBatchValidationLike>} | undefined;
  // once per source and mapping: the importer can be slow, and the preview and the commit read the same plan
  private _planned: {frame: DG.DataFrame, key: string, plan: Promise<Planned[]>} | undefined;
  // an older dry run answering last must not paint over a newer one
  private _previewGen = 0;
  // '' while the dry run runs, a reason where all-or-nothing would abort
  private readonly _previewBlock = signal<string | null>(null);
  private readonly _importText = signal('Import');
  private readonly _aborted = signal(false);

  constructor(private readonly _table: DomainTable, private readonly _targets: IProperty[],
    private readonly _options: DomainImportOptions,
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
    if (batch.upsert) {
      this._mode = new ChoiceInput({label: 'Mode', name: 'mode', value: 'insert', nullable: false,
        items: [{value: 'insert', label: 'Add new rows'}, {value: 'upsert', label: 'Add or update'}]});
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
    return divV([table.root, statusLine(signal(null), computed(() => this._mappingStatus()))],
      'u2-domain-import-mapping-step');
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
      .map((name) => `${ImportFlow._caption(this, name)} ${mapped.has(name) ? '✓' : '✗'}`);
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
    try {
      plan = await this._plan(frame, mapping,
        (done) => status.textContent = `Checking rows: ${done} of ${frame.rowCount}`);
      const sent = ImportFlow._sent(plan);
      if (this._table.table.validate !== undefined) {
        const key = JSON.stringify([mapping, options]);
        if (this._checked?.key !== key) {
          const scanned = sent.slice(0, PREVIEW_ROWS).map((i) => plan[i].payload!);
          this._checked = {key, report: this._table.table.validate(scanned, options)};
        }
        report = await this._checked.report;
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
    if (report === null)
      this._previewBlock.value = null;
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
    const kept = plan.filter((p) => p.payload === null);
    const keptErrors = kept.filter((p) => p.verdict.error !== undefined).length;
    const errors = report.errorCount + keptErrors;
    const checked = sent.length > PREVIEW_ROWS ? `Checked the first ${PREVIEW_ROWS.toLocaleString()} of ` +
      `${sent.length.toLocaleString()} rows — the rest are checked on import:` : '';
    host.prepend(ImportFlow._headline(ImportFlow._counts([[report.willInsert, 'to add', 'success'],
      [report.willUpdate, 'to update', 'success'], [report.willSkip + kept.length - keptErrors, 'to skip', 'default'],
      [errors, errors === 1 ? 'error' : 'errors', 'error']]), 'Nothing to import'));
    if (checked !== '')
      host.prepend(span(checked, 'u2-domain-import-summary'));
    if (showing !== '')
      host.append(span(`${showing}.`, 'u2-domain-import-summary'));
    const valid = sent.length > PREVIEW_ROWS ? null : report.willInsert + report.willUpdate;
    const noun = (n: number) => plural(n, nounOf(this._table.info, true), nounOf(this._table.info));
    this._importText.value = `Import ${noun(valid ?? sent.length)}`;
    // the counts are what the batch WOULD do row by row; all-or-nothing makes one bad row the
    // verdict on all of them
    this._previewBlock.value = options.allOrNothing !== false && errors > 0 ? this._abortReason(errors, valid) : null;
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
    return {index, values: planned.payload ?? planned.cells, kind, label: planned.verdict.label ?? '', note, errors,
      id: null, key: ''};
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
    return [ImportFlow._numberColumn(), ImportFlow._resultColumn(VERDICTS),
      ...mapping.map(([target]) => this._valueColumn(target)), ImportFlow._noteColumn()];
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
    return target === undefined ? '' : link(row.key, () => {
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
    const rows = this._payload(frame, mapping, frame.rowCount);
    const prepare = this._table.importer?.prepare;
    const verdicts = prepare === undefined ? [] : await prepare(rows, progress);
    return rows.map((cells, i): Planned => {
      const verdict = verdicts[i] ?? {};
      if (verdict.error !== undefined || verdict.existing !== undefined || verdict.sameAs !== undefined)
        return {cells, payload: null, verdict};
      const payload = {...cells, ...verdict.values};
      const problem = this._validate(payload);
      return problem === null ? {cells, payload, verdict} : {cells, payload: null, verdict: {error: problem}};
    });
  }

  private _validate(row: Record<string, unknown>): string | null {
    const validators = this._table.validators;
    for (const column of validators.columns) {
      const problem = validators.check(column, row[column], row as RowView);
      if (problem !== null)
        return `${ImportFlow._caption(this, column)}: ${problem}`;
    }
    return null;
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
      for (let i = 0; i < count; i++)
        rows[i][target] = ImportFlow._valueOrNull(column, i);
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
    for (const choice of this._choices.values())
      choice.dispose();
    this._choices.clear();
    if (frame === null) {
      this._rows.value = [];
      return;
    }
    // the asterisk the forms use, in the option itself: a picker's own list is all the user sees
    const items: ChoiceItem[] = [{value: SKIP, label: SKIP_LABEL}, ...this._targets.map((p) =>
      ({value: p.name!, label: `${DomainTable.captionOf(p)}${p.nullable === false ? ' *' : ''}`}))];
    const columns = frame.columns.names();
    for (const column of columns) {
      // owned by hand, not by a scope: the first remap runs inside the wizard's own constructor
      // (a source handed in as an option), where there is no wizard to own anything yet
      const choice = new ChoiceInput({items, value: this._match(frame.columns.byName(column)), inline: true,
        name: column, nullable: false});
      // required as the forms say it: a target the table will not take empty
      choice.effect(() => {
        const target = choice.value;
        choice.required = this._targets.some((p) => p.name === target && p.nullable === false);
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
    const prop = flow._table.properties.find((p) => p.name === column);
    return prop === undefined ? column : DomainTable.captionOf(prop);
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
    return hit((a) => a === want) ?? hit((a) => a.startsWith(want) || want.startsWith(a)) ??
      hit((a) => a.includes(want) || want.includes(a)) ?? (typed.length === 1 ? typed[0].name! : SKIP);
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
    const upsert = this._mode?.peek() === 'upsert';
    for (const name of this._required()) {
      if (seen.has(name))
        continue;
      const caption = ImportFlow._caption(this, name);
      problems.push(upsert ? `Add or update matches rows by ${caption} — map a column to "${caption}".` :
        `Required column "${caption}" is not mapped — rows without a value will fail.`);
    }
    return problems;
  }

  private _required(): string[] {
    return this._mode?.peek() === 'upsert' ? this._table.info.businessKey : this._targets
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
      const report = refused ? {inserted: 0, updated: 0, skipped: 0, errorCount: 0, error: 'validation', rows: []} :
        await this._table.table.batch!(sent.map((i) => plan[i].payload!), options);
      this._result = ImportFlow._merged(report, sent, plan);
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

  private static _merged(report: DomainBatchReportLike, sent: number[], plan: Planned[]): DomainBatchReportLike {
    const kept = plan.flatMap((p, i): DomainBatchReportLike['rows'] => p.payload !== null ? [] :
      [p.verdict.error !== undefined ? {index: i, id: null, status: 'error', errors: [{message: p.verdict.error}]} :
        {index: i, id: null, status: 'skipped', existingId: p.verdict.existing}]);
    const errors = kept.filter((r) => r.status === 'error').length;
    return {...report, skipped: (report.skipped ?? 0) + kept.length - errors,
      errorCount: report.errorCount === undefined ? undefined : report.errorCount + errors,
      rows: [...report.rows.map((r) => ({...r, index: sent[r.index]})), ...kept].sort((a, b) => a.index - b.index)};
  }

  private async _render(report: DomainBatchReportLike, plan: Planned[], frame: DG.DataFrame): Promise<void> {
    const failed = report.error !== undefined;
    // a warehouse refusal carries the failing rows and no totals
    const errorCount = report.errorCount ?? report.rows.filter((row) => row.errors?.length).length;
    const counts = ImportFlow._counts(ImportFlow._landed(report, errorCount));
    const aborted = `Import aborted — ${plural(errorCount, 'row has errors', 'rows have errors')}; ` +
      'nothing was committed.';
    this._report.replaceChildren(ImportFlow._headline(failed ? [[aborted, 'error']] : counts, 'Nothing was imported'));
    this._report.classList.toggle('u2-domain-import-aborted', failed);
    this._wizard.setStepFailed('report', failed || errorCount > 0);
    const outcomes = this._outcomes(report, plan, failed);
    let saved = new Map<string, RowView>();
    try {
      saved = await this._rowsById([...new Set(outcomes.map((o) => o.id).filter((id): id is string => id !== null))]);
    } catch (e) {
      this._report.append(span(`The keys could not be read: ${DomainErrors.message(e)}`, 'u2-domain-import-problem'));
    }
    for (const row of outcomes)
      row.key = row.id === null || !saved.has(row.id) ? '' : this._keyOf(saved.get(row.id)!);
    // a dialog closed while the keys were read has no report to show; the write-back and the balloon still follow
    if (!this._wizard.scope.isDisposed)
      this._report.append(this._rowsTable(this._wizard, 'report', outcomes, this._reportColumns(saved)));
    const column = this._keyColumn();
    if (!failed && column !== null && this._writeBack?.peek()) {
      if (frame.columns.names().includes(column))
        frame.columns.remove(column);
      // column names ignore case: a source's own "sku" keeps its values beside the "Sku" written back
      const added = frame.columns.addNewString(frame.columns.getUnusedName(column)).init((i) => outcomes[i].key);
      this._report.append(span(`Added a "${added.name}" column to ${frame.name} with each row's ${this._keyHeader()}.`,
        'u2-domain-import-summary'));
    }
    if (failed) {
      const refused = report.rows.filter((r) => r.status === 'error' && plan[r.index].payload !== null).length;
      const valid = ImportFlow._sent(plan).length - refused;
      if (this._allOrNothing?.peek() && valid > 0) {
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

  // a repeat takes its first row's id
  private _outcomes(report: DomainBatchReportLike, plan: Planned[], failed: boolean): DomainImportRow[] {
    const bySource = new Map(report.rows.map((r) => [r.index, r]));
    const outcomes = plan.map((p, i): DomainImportRow => {
      const row = bySource.get(i);
      const errors = row?.errors ?? [];
      if (row?.status === 'error')
        return this._row(i, p, 'error', p.payload === null ? p.verdict.error ?? '' : this._describe(errors), errors);
      if (row?.status === 'duplicate')
        return {...this._row(i, p, 'skip', this._duplicate(), errors), id: row.existingId ?? null};
      if (p.payload === null)
        return this._kept(i, p);
      if (failed || row === undefined)
        return this._row(i, p, failed ? 'aborted' : '', failed ? 'Another row has errors' : '', errors);
      return {...this._row(i, p, row.status === 'inserted' ? 'insert' : 'update', p.verdict.note ?? '', errors),
        id: row.id};
    });
    plan.forEach((p, i) => {
      if (p.verdict.sameAs !== undefined)
        outcomes[i].id = outcomes[p.verdict.sameAs]?.id ?? null;
    });
    return outcomes;
  }

  private async _rowsById(ids: string[]): Promise<Map<string, RowView>> {
    const rows = new Map<string, RowView>();
    for (let at = 0; at < ids.length; at += KEY_CHUNK) {
      const chunk = ids.slice(at, at + KEY_CHUNK);
      const filter = {property: 'id', operator: '=', value: chunk};
      for (const row of await this._table.table.query({filter, limit: chunk.length}))
        rows.set(String(row.id), row as RowView);
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

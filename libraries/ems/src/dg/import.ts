/* `domains.import` — the import wizard (WO 3-9), the u2 half of `domain_import_dialog.dart`:
   **source** (any open frame, or any file the platform can read, through `ui.input.table`),
   **mapping** (a `DataTable` of source column → target, auto-matched by name and semantic type; only
   the columns the caller may write are offered), **preview** (the server's dry run over the first 1000
   rows — `validate`, which is the commit's own checks inside a rolled-back transaction; the mapped
   rows alone where the storage has no dry run) and **report** (every row's outcome and key). What is
   posted is the MAPPED columns only, under their target names, plus what the table's
   {@link DomainImporter} derives, so a skipped or renamed source column never reaches the server;
   the source frame is touched only to write the keys back, when asked. The options offered are
   the ones `support.batch` declares: an option the storage refuses is not on the form and not sent.
   An app makes the import its own through `table.importer`: the rows it derives or keeps out, the
   columns of the preview and the report, how a refusal reads. */
import * as ui from 'datagrok-api/ui';
import type * as DG from 'datagrok-api/dg';
import {signal, Control, button, div, divV, span, plural, text, link, Wizard, DataTable, ChoiceInput,
  BoolInput, badge, Form, notify} from '@datagrok-libraries/u2';
import type {IProperty, ChoiceItem, BadgeVariant, RowView} from '@datagrok-libraries/u2';
import type {DomainBatchOptionsLike, DomainBatchReportLike, DomainBatchValidationLike} from '../core/domain-backend.js';
import {fromDartInput} from '@datagrok-libraries/u2/src/dg/index.js';
import {DomainTable} from './facade.js';
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
 * index) or as an `error`. `label` names it on the row's badge, `note` says why. */
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
  /** The row it is now — in the report, once read back — by id and as the table captions it. */
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
/** The preview scans at most this many rows — it is advisory, and the server re-validates every
 * row it is sent; the commit itself is unbounded. */
const PREVIEW_ROWS = 1000;
/** Rows the preview draws — what the import looks like, not all of it. */
const PREVIEW_SHOWN = 20;
/** A row's kind as the preview shows it: the verb, not the wire's op name. */
const VERDICTS: Record<string, {label: string, variant: BadgeVariant}> = {
  insert: {label: 'Add', variant: 'accent'},
  update: {label: 'Update', variant: 'accent'},
  skip: {label: 'Skip', variant: 'default'},
  error: {label: 'Error', variant: 'error'},
};
/** A row's kind as the report shows it, once the batch has run. */
const OUTCOMES: Record<string, {label: string, variant: BadgeVariant}> = {
  insert: {label: 'Added', variant: 'accent'},
  update: {label: 'Updated', variant: 'accent'},
  skip: {label: 'Skipped', variant: 'default'},
  error: {label: 'Error', variant: 'error'},
  aborted: {label: 'Not imported', variant: 'default'},
};
/** Why a row whose business key the table or the file already has is skipped. */
const DUPLICATE = 'Duplicate business key';
/** Rows of the report's key lookup per query. */
const KEY_CHUNK = 500;

/** One row of the mapping step — the target is the choice input the flow keeps for that column. */
interface Mapping {column: string}

/** Opens the wizard over `table`; resolves to the server's report, or to null where it was
 * cancelled. Exported as `domains.import` — `import` cannot name a function. */
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
  /** Bumped by every edit the gates read — the source, the mapping and the mode. */
  private readonly _recheck = signal(0);
  private readonly _report = div([], 'u2-domain-import-report');
  /** Absent where `support.batch` refuses the option. */
  private _mode?: ChoiceInput;
  private _allOrNothing?: BoolInput;
  private _errorOnDuplicate?: BoolInput;
  /** Absent where no source column takes the keys back. */
  private _writeBack?: BoolInput;
  private _wizard!: Wizard;
  private _committed = false;
  private _result: DomainBatchReportLike | null = null;
  private _settled = false;
  /** The dry run and what it was run for — the mapping, the mode and the flags; the preview posts
   * again only when that key changed. */
  private _checked: {key: string, report: Promise<DomainBatchValidationLike>} | undefined;
  /** Every row's verdict and what it sends, made once per source and mapping: the importer can be
   * slow, and the preview and the commit read the same plan. */
  private _planned: {frame: DG.DataFrame, key: string, plan: Promise<Planned[]>} | undefined;
  /** Bumped by every post: an older dry run answering last must not paint over a newer one. */
  private _previewGen = 0;
  /** What the dry run says of committing: '' blocks while it runs, a reason where all-or-nothing
   * would abort. */
  private readonly _previewBlock = signal<string | null>(null);
  private readonly _importText = signal('Import');

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
          onActivate: () => void this._commit()},
      ],
      onFinish: () => this._settle(this._result),
      onCancel: () => this._settle(null),
    });
    this._wizard.root.classList.add('u2-domain-import');
    this._wizard.openInDialog(`Import ${ImportFlow._plural(this._table)}`, {width: 720, height: 560});
  }

  /** The first blocking problem, which is what the NEXT button shows and gates on. */
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
    if (this._options.source !== undefined)
      picker.value.value = this._options.source;
    const summary = span('', 'u2-domain-import-summary');
    const batch = this._table.table.support.batch;
    if (batch.upsert) {
      this._mode = new ChoiceInput({label: 'Mode', name: 'mode', value: 'insert', nullable: false,
        items: [{value: 'insert', label: 'Add new rows'}, {value: 'upsert', label: 'Add or update'}]});
    }
    if (batch.partial) {
      this._allOrNothing = new BoolInput({label: 'All or nothing', value: true, name: 'allOrNothing',
        tooltipText: 'Any invalid row aborts the whole import'});
    }
    if (batch.skipDuplicates) {
      this._errorOnDuplicate = new BoolInput({label: 'Error on duplicate', value: false, name: 'errorOnDuplicate',
        tooltipText: 'Treat business-key duplicates as errors instead of skips (insert mode)'});
    }
    const column = this._keyColumn();
    if (column !== null) {
      this._writeBack = new BoolInput({label: `Add a "${column}" column`, value: true, name: 'writeBack',
        tooltipText: `Each row of the source table gets the ${this._keyHeader()} it was imported as, or found as`});
    }
    const form = new Form().addAll([picker, this._mode, this._allOrNothing, this._errorOnDuplicate, this._writeBack]
      .filter((input) => input !== undefined));
    form.effect(() => {
      const frame = picker.value.value;
      summary.textContent = frame === null ? '' :
        `${frame.name}: ${frame.rowCount} rows, ${frame.columns.names().length} columns`;
      this._frame.value = frame;
      this._remap(frame);
      this._touch();
    });
    form.effect(() => {
      const insert = this._mode?.value.value !== 'upsert';
      if (this._errorOnDuplicate !== undefined)
        this._errorOnDuplicate.enabled = insert;
      this._touch();
    });
    return divV([form.root, summary], 'u2-domain-import-source');
  }

  private _mappingStep(): HTMLElement {
    const table = new DataTable<Mapping>({
      columns: [
        {name: 'column', header: 'Source column'},
        {name: 'target', header: 'Target column', render: (item) => this._choices.get(item.column)!.root},
      ],
      keyOf: (item) => item.column,
      rowHeight: 32,
      items: this._rows,
    });
    table.root.classList.add('u2-domain-import-mapping');
    return table.root;
  }

  /** The table importer's verdicts over every row, then the server's own over the first
   * {@link PREVIEW_ROWS} rows it would be sent: the dry run runs every check the commit runs,
   * inside a transaction it rolls back, so the preview and the report cannot disagree. Posted once
   * per entry into the step, and again only when the mapping, the mode or the source changed since. */
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
      // the commit is the authority: a preview that cannot reach the server says so and lets the
      // import through
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
    if (report === null) {
      const lead = rows.length < frame.rowCount ? `Showing the first ${rows.length.toLocaleString()} of ` +
        `${frame.rowCount.toLocaleString()} rows — rows` : 'Rows';
      host.append(span(`${lead} are checked when imported.`, 'u2-domain-import-summary'));
      return;
    }
    const kept = plan.filter((p) => p.payload === null);
    const keptErrors = kept.filter((p) => p.verdict.error !== undefined).length;
    const errors = report.errorCount + keptErrors;
    host.append(span(`${report.willInsert} will be added, ${report.willUpdate} updated, ` +
      `${report.willSkip + kept.length - keptErrors} skipped, ${plural(errors, 'row has errors', 'rows have errors')}.`,
    'u2-domain-import-summary'));
    const valid = sent.length > PREVIEW_ROWS ? null : report.willInsert + report.willUpdate;
    this._importText.value = valid === null ? 'Import' :
      `Import ${plural(valid, ImportFlow._singular(this._table), ImportFlow._plural(this._table))}`;
    // the counts are what the batch WOULD do row by row; all-or-nothing makes one bad row the
    // verdict on all of them
    this._previewBlock.value = options.allOrNothing !== false && errors > 0 ? this._abortReason(errors, valid) : null;
    if (sent.length > PREVIEW_ROWS) {
      host.append(span(`Checked the first ${PREVIEW_ROWS.toLocaleString()} of ` +
        `${sent.length.toLocaleString()} rows — the rest are checked on import.`,
      'u2-domain-import-summary'));
    }
  }

  /** Why the preview will not import: all-or-nothing would abort over the rows with errors. `valid`
   * is null where the dry run did not see every row. */
  private _abortReason(errors: number, valid: number | null): string {
    const lead = plural(errors, 'row has errors', 'rows have errors');
    if (this._allOrNothing === undefined)
      return `${lead}, and this import is all or nothing.`;
    if (valid === 0)
      return `${lead}, and no row is valid.`;
    return `${lead}. Uncheck "All or nothing" on Source to import the ` +
      `${valid === null ? 'valid rows' : plural(valid, 'valid row', 'valid rows')} and skip the rest.`;
  }

  /** Every source row as the preview shows it: the importer's verdict where it keeps the row out,
   * the dry run's where the row was sent to it. */
  private _previewRows(plan: Planned[], sent: number[], report: DomainBatchValidationLike | null): DomainImportRow[] {
    const rows = plan.map((p, i) => p.payload === null ? this._kept(i, p) :
      this._row(i, p, report === null ? 'insert' : '', p.verdict.note ?? '', []));
    for (const dry of report?.rows ?? []) {
      const at = sent[dry.index];
      const errors = dry.errors ?? [];
      rows[at] = this._row(at, plan[at], dry.predicted, errors.length > 0 ? this._describe(errors) :
        plan[at].verdict.note ?? (dry.predicted === 'skip' ? DUPLICATE : ''), errors);
    }
    return rows;
  }

  /** A row the plan keeps out of the batch. */
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

  /** The rows as a table, its columns the importer's choice; the cells a refusal names are marked. */
  private _rowsTable(owner: Control, step: DomainImportStep, rows: DomainImportRow[],
    defaults: DomainImportColumn[]): HTMLElement {
    const columns = this._table.importer?.columns?.(defaults, step) ?? defaults;
    const byIndex = new Map(rows.map((row) => [row.index, row]));
    const table = owner.runInScope(() => new DataTable<DomainImportRow>({
      rowHeight: this._table.importer?.rowHeight ?? 24,
      keyOf: (row) => String(row.index),
      columns: columns.map((c) => ({name: c.name, header: c.header, width: c.width,
        render: (row, _index, cell) => ImportFlow._cell(c.render(row), cell)})),
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

  /** A cell's content, its text the cell's title, since a narrow column cuts it; plain text goes
   * in a span of its own, which can trail off inside the cell that fills the row. */
  private static _cell(content: HTMLElement | string, cell: HTMLElement): HTMLElement {
    const label = typeof content === 'string' ? content : content.textContent ?? '';
    if (label !== '')
      cell.title = label;
    return typeof content === 'string' ? span(content, 'u2-domain-import-text') : content;
  }

  /** What the preview shows by default: the row's number and verdict, its mapped cells, the note. */
  private _previewColumns(mapping: [string, string][]): DomainImportColumn[] {
    return [ImportFlow._numberColumn('#', '44px'), ImportFlow._resultColumn(VERDICTS),
      ...mapping.map(([target]) => this._valueColumn(target)),
      {name: '~note', header: 'Note', render: (row) => row.note}];
  }

  /** What the report shows by default: the row's number and outcome, the row it is now, why. */
  private _reportColumns(saved: Map<string, RowView>): DomainImportColumn[] {
    return [ImportFlow._numberColumn('Source row', '84px'), ImportFlow._resultColumn(OUTCOMES),
      {name: '~key', header: this._keyHeader(), width: '120px', render: (row) => this._keyLink(row, saved)},
      {name: '~note', header: 'Message', render: (row) => row.note}];
  }

  private _valueColumn(target: string): DomainImportColumn {
    return {name: target, header: ImportFlow._caption(this, target), render: (row) => text(row.values[target])};
  }

  private static _numberColumn(header: string, width: string): DomainImportColumn {
    return {name: '~row', header, width, render: (row) => String(row.index + 1)};
  }

  private static _resultColumn(kinds: Record<string, {label: string, variant: BadgeVariant}>): DomainImportColumn {
    return {name: '~result', header: 'Result', width: '110px', render: (row) => {
      const kind = kinds[row.kind];
      const label = row.label || kind?.label;
      return label === undefined ? '' : badge(label, {variant: kind?.variant ?? 'default'});
    }};
  }

  /** The row an outcome is about, as a link: the wizard closes and the row opens. */
  private _keyLink(row: DomainImportRow, saved: Map<string, RowView>): HTMLElement | string {
    const target = row.id === null ? undefined : saved.get(row.id);
    return target === undefined ? '' : link(row.key, () => {
      this._settle(this._result);
      this._table.open(target);
    });
  }

  /** Refusals as a row's note, each named by the field it is about. */
  private _describe(errors: DomainImportError[]): string {
    return errors.map((e) => {
      const caption = e.column ? ImportFlow._caption(this, e.column) : '';
      return this._table.importer?.describe?.(e, caption) ?? (e.code === 'unique' && caption ?
        `${caption} is already used` : caption ? `${caption}: ${e.message}` : e.message);
    }).join('; ');
  }

  /** Every row's verdict and payload — the importer's verdicts, then the table's validators over
   * what would be sent: an import is held to the rules a form holds a row to. */
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

  /** The source rows the batch is sent, in order — a sent row's place in the batch is its index here. */
  private static _sent(plan: Planned[]): number[] {
    return plan.flatMap((p, i) => p.payload === null ? [] : [i]);
  }

  /** The mapped rows, under their target names — what is posted, and the only thing that is. */
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

  /** Only the options the form offers; an absent one is left to the server's default. */
  private _batchOptions(): DomainBatchOptionsLike {
    const mode = this._mode?.value.peek() === 'upsert' ? 'upsert' : 'insert';
    const options: DomainBatchOptionsLike = {mode};
    if (this._allOrNothing !== undefined)
      options.allOrNothing = this._allOrNothing.value.peek();
    if (this._errorOnDuplicate !== undefined)
      options.errorOnDuplicate = mode === 'insert' && this._errorOnDuplicate.value.peek();
    return options;
  }

  private _touch(): void {
    this._recheck.value = this._recheck.peek() + 1;
  }

  /** One target choice per source column — the importer's match, else the auto-match. The choices
   * are kept here, not in the mapping table's cells: the table recycles its rows, and a choice that
   * lived in a cell would be rebuilt — and lose its value — on every scroll. */
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
    const items: ChoiceItem[] = [SKIP, ...this._targets.map((p) =>
      ({value: p.name!, label: `${ImportFlow.captionOf(p)}${p.nullable === false ? ' *' : ''}`}))];
    const columns = frame.columns.names();
    for (const column of columns) {
      // owned by hand, not by a scope: the first remap runs inside the wizard's own constructor
      // (a source handed in as an option), where there is no wizard to own anything yet
      const choice = new ChoiceInput({items, value: this._match(frame.columns.byName(column)), inline: true,
        name: column, nullable: false});
      // the asterisk the forms use: a target the table will not take empty
      choice.effect(() => choice.root.classList.toggle('u2-input-required',
        this._targets.some((p) => p.name === choice.value.value && p.nullable === false)));
      choice.value.subscribe(() => this._touch());
      this._choices.set(column, choice);
    }
    this._rows.value = columns.map((column) => ({column}));
  }

  /** A target as a form labels it: the declared caption, else the column name made presentable. */
  static captionOf(prop: IProperty): string {
    const name = prop.friendlyName ?? (prop.name ?? '').replace(/_/g, ' ');
    return `${name.charAt(0).toUpperCase()}${name.slice(1)}`;
  }

  /** The same, by column name. */
  private static _caption(flow: ImportFlow, column: string): string {
    const prop = flow._table.properties.find((p) => p.name === column);
    return prop === undefined ? column : ImportFlow.captionOf(prop);
  }

  /** The source column the keys are written back to: the importer's choice, else the business
   * key's caption where the key is one column. */
  private _keyColumn(): string | null {
    const chosen = this._table.importer?.keyColumn;
    const key = this._table.info.businessKey;
    return chosen === false ? null : chosen ?? (key.length === 1 ? ImportFlow._caption(this, key[0]) : null);
  }

  /** What the report's key column is headed: the business key's caption, else the table's noun. */
  private _keyHeader(): string {
    const key = this._table.info.businessKey;
    if (key.length === 1)
      return ImportFlow._caption(this, key[0]);
    const name = (this._table.info.singularName || 'row').replace(/_/g, ' ');
    return `${name.charAt(0).toUpperCase()}${name.slice(1)}`;
  }

  /** The table's singular name, as the user reads it. */
  private static _singular(table: DomainTable): string {
    return (table.info.singularName || 'row').replace(/_/g, ' ').toLowerCase();
  }

  /** The table's plural name, as the user reads it. */
  private static _plural(table: DomainTable): string {
    return (table.info.pluralName || 'rows').replace(/_/g, ' ').toLowerCase();
  }

  /** Where a source column goes: the importer's match, else the name the way a person reads it —
   * "Chemical name" is the name column and "CAS number" the cas one, so the comparison drops case and
   * everything that is not a letter or a digit, and an exact hit anywhere beats a partial hit
   * somewhere earlier — else the one target of the column's semantic type. */
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

  /** Target column → source column, in mapping order; a target claimed twice keeps the first. */
  private _mapping(): Map<string, string> {
    const out = new Map<string, string>();
    for (const [column, choice] of this._choices) {
      const target = choice.value.peek();
      if (target !== null && target !== SKIP && !out.has(target))
        out.set(target, column);
    }
    return out;
  }

  /** Everything that blocks the import, in the order `domain_import_dialog.dart` reports it. */
  private _problems(): string[] {
    this._recheck.value;
    if (this._frame.value === null)
      return ['Choose an open table or a file'];
    const problems: string[] = [];
    const seen = new Set<string>();
    for (const choice of this._choices.values()) {
      const target = choice.value.peek();
      if (target === null || target === SKIP)
        continue;
      if (seen.has(target))
        problems.push(`Two source columns are mapped to "${target}" — keep one and skip the other.`);
      seen.add(target);
    }
    if (seen.size === 0)
      problems.push('Map at least one column.');
    if (this._mode?.value.peek() === 'upsert') {
      for (const key of this._table.info.businessKey) {
        if (!seen.has(key))
          problems.push(`Upsert merges by the business key — map a column to "${ImportFlow._caption(this, key)}".`);
      }
    } else {
      for (const p of this._targets) {
        if (p.nullable === false && (p.defaultValue === undefined || p.defaultValue === null) && !seen.has(p.name!)) {
          problems.push(`Required column "${ImportFlow.captionOf(p)}" is not mapped — ` +
            'rows without a value will fail.');
        }
      }
    }
    return problems;
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
      notify.error(DomainErrors.message(e));
      return;
    }
    await this._render(this._result, plan, frame);
  }

  /** The batch's report over the source rows: its rows under their source numbers, and the rows
   * the plan kept out of it as skips and errors of their own — what the caller gets back. */
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
    this._report.replaceChildren(span(failed ?
      `Import aborted — ${plural(errorCount, 'row has errors', 'rows have errors')}; ` +
      'nothing was committed.' :
      `${ImportFlow._landed(report)}, ${report.skipped ?? 0} skipped, ${errorCount} failed.`,
    failed ? 'u2-domain-import-problem' : 'u2-domain-import-summary'));
    const outcomes = this._outcomes(report, plan, failed);
    let saved = new Map<string, RowView>();
    try {
      saved = await this._rowsById([...new Set(outcomes.map((o) => o.id).filter((id): id is string => id !== null))]);
    } catch (e) {
      this._report.append(span(`The keys could not be read: ${DomainErrors.message(e)}`, 'u2-domain-import-problem'));
    }
    for (const row of outcomes)
      row.key = row.id === null || !saved.has(row.id) ? '' : this._table.renderer.caption(saved.get(row.id)!);
    this._report.append(this._rowsTable(this._wizard, 'report', outcomes, this._reportColumns(saved)));
    const column = this._keyColumn();
    if (!failed && column !== null && this._writeBack?.value.peek()) {
      if (frame.columns.names().includes(column))
        frame.columns.remove(column);
      frame.columns.addNewString(column).init((i) => outcomes[i].key);
      this._report.append(span(`Each row's ${this._keyHeader()} is in the "${column}" column of ${frame.name}.`,
        'u2-domain-import-summary'));
    }
    if (failed) {
      const refused = report.rows.filter((r) => r.status === 'error' && plan[r.index].payload !== null).length;
      this._report.append(...this._wayBack(ImportFlow._sent(plan).length - refused));
      notify.error('Import failed — nothing was committed. See the report.');
    }
    else
      notify.info(`Imported: ${ImportFlow._landed(report)}`);
  }

  /** Every source row's outcome: the batch's, else the plan's; a repeat takes its first row's id. */
  private _outcomes(report: DomainBatchReportLike, plan: Planned[], failed: boolean): DomainImportRow[] {
    const bySource = new Map(report.rows.map((r) => [r.index, r]));
    const outcomes = plan.map((p, i): DomainImportRow => {
      const row = bySource.get(i);
      const errors = row?.errors ?? [];
      if (row?.status === 'error')
        return this._row(i, p, 'error', p.payload === null ? p.verdict.error ?? '' : this._describe(errors), errors);
      if (row?.status === 'duplicate')
        return {...this._row(i, p, 'skip', DUPLICATE, errors), id: row.existingId ?? null};
      if (p.payload === null)
        return this._kept(i, p);
      if (failed || row === undefined)
        return this._row(i, p, failed ? 'aborted' : '', '', errors);
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

  /** What landed; `merged` stands in for `updated` where the storage cannot tell the two apart. */
  private static _landed(report: DomainBatchReportLike): string {
    const merged = report.merged ?? 0;
    const updated = report.updated ?? 0;
    const counts = [`${report.inserted ?? 0} inserted`];
    if (merged === 0 || updated > 0)
      counts.push(`${updated} updated`);
    if (merged > 0)
      counts.push(`${merged} merged`);
    return counts.join(', ');
  }

  /** The way out of an aborted all-or-nothing run: what to change, and the BACK the report step
   * has none of — the mapping and the source stand, so only the flag has to be turned off. */
  private _wayBack(valid: number): HTMLElement[] {
    const out: HTMLElement[] = [];
    if (this._allOrNothing?.value.peek() && valid > 0) {
      out.push(span(`Uncheck "All or nothing" to import the ${plural(valid, 'valid row', 'valid rows')} ` +
        'and skip the rest.', 'u2-domain-import-summary'));
    }
    out.push(div([button('BACK', () => {
      this._committed = false;
      this._result = null;
      this._checked = undefined;
      this._wizard.goTo('source');
    })], 'u2-domain-import-actions'));
    return out;
  }
}

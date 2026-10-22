/* `domains.import` — the import wizard (WO 3-9), the u2 half of `domain_import_dialog.dart`:
   **source** (any open frame, or any file the platform can read, through `ui.input.table`),
   **mapping** (a `DataTable` of source column → target, auto-matched by name, else by semType; only
   the columns the caller may write are offered), **preview** (the server's dry run over the first
   1000 rows — `validate`, which is the commit's own checks inside a rolled-back transaction; the
   mapped rows alone where the storage has no dry run) and **report** (the commit's, each row that
   landed by its key). What is posted is the MAPPED columns only, under their target names — a
   reference column's names as the ids of the rows they name — so a skipped or renamed source column
   never reaches the server and the source frame is never touched; a row the import found a problem
   with, or a derivation skipped (`DomainSkip`), is not posted at all. `domains.importRows` is the
   same run without the wizard. The options offered are the ones `support.batch` declares: an option
   the storage refuses is not on the form and not sent. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import type * as DG from 'datagrok-api/dg';
import {signal, Control, button, div, divV, span, link, plural, text, isEmpty, Filters, Wizard, DataTable, BasicTable,
  ChoiceInput, BoolInput, badge, Form, notify} from '@datagrok-libraries/u2';
import type {IProperty, ChoiceItem, BadgeVariant, RowView, FilterGroup} from '@datagrok-libraries/u2';
import type {DomainBatchOptionsLike, DomainBatchReportLike, DomainBatchValidationLike,
  DomainBatchValidationRowLike} from '../core/domain-backend.js';
import {fromDartInput, moleculeRenderer} from '@datagrok-libraries/u2/src/dg/index.js';
import {DomainSource} from '../core/domain-source.js';
import {DomainAddress} from './address.js';
import {DomainTable} from './facade.js';
import type {RowProblem} from './facade.js';
import {DomainErrors} from './errors.js';
import {DomainPick} from './pick.js';

export interface DomainImportOptions {
  /** The frame to import, when the caller already has one; the picker still offers the rest. */
  source?: DG.DataFrame;
  /** What "Show the N imported rows" does, once the wizard has closed, with the query naming them. */
  show?(query: FilterGroup): void;
}

/** The columns a report frame adds to the source rows; the key column is named as the table's caption column. */
export const IMPORT_COLUMNS = {result: 'Import result', message: 'Import message'} as const;

/** The preview scans at most this many rows — it is advisory, and the server re-validates every
 * row it is sent; the commit itself is unbounded. */
const PREVIEW_ROWS = 1000;
/** Issue lines rendered in a preview or a report; the counts stay exact. */
const ISSUE_CAP = 20;
/** Rows the preview draws — what the import looks like, not all of it. */
const PREVIEW_SHOWN = 20;
/** The preview's leading column: what the server said it would do with the row. */
const VERDICT = '~predicted';
/** A verdict as the preview shows it: the verb, not the wire's op name. */
const VERDICTS: Record<string, {label: string, variant: BadgeVariant}> = {
  insert: {label: 'Add', variant: 'accent'},
  update: {label: 'Update', variant: 'accent'},
  skip: {label: 'Skip', variant: 'default'},
  error: {label: 'Error', variant: 'error'},
};
/** What became of a row, by the status its report line has. */
const OUTCOMES: Record<string, string> = {inserted: 'Added', updated: 'Updated', merged: 'Merged',
  duplicate: 'Skipped', skipped: 'Skipped', error: 'Error'};
const LANDED = ['inserted', 'updated', 'merged'];
/** What a source column maps to when it maps to nothing. */
const SKIP = '(skip)';
const THUMBNAIL = moleculeRenderer({width: 80, height: 40});

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
  const targets = ImportRun.targets(table);
  if (targets.length === 0) {
    notify.warning(`There is no column of ${table.address} you may write.`);
    return Promise.resolve(null);
  }
  return new Promise((resolve) => new ImportFlow(table, targets, options, resolve).open());
}

/** The import without the wizard: the columns matched as the wizard matches them, `mapping` (target → source column)
 * over that; resolves to the frame "Open as table" opens. */
export async function importRows(table: DomainTable, frame: DG.DataFrame,
  options: DomainBatchOptionsLike & {mapping?: Record<string, string>} = {}): Promise<DG.DataFrame> {
  const {mapping = {}, ...batch} = options;
  const targets = ImportRun.targets(table);
  const pairs = new Map<string, string>();
  for (const column of frame.columns.toList()) {
    const target = ImportRun.autoMatch(targets, column);
    if (target !== undefined && !pairs.has(target) && !Object.values(mapping).includes(column.name))
      pairs.set(target, column.name);
  }
  for (const [target, column] of Object.entries(mapping))
    pairs.set(target, column);
  const run = new ImportRun(table, frame, [...pairs]);
  await run.commit(batch);
  return run.reportFrame();
}

/** One import of a frame's mapped columns: the wizard's preview and commit, and {@link importRows}. */
class ImportRun {
  readonly rows: Record<string, unknown>[];
  problems: RowProblem[] = [];
  /** The source row of each row sent: the batch numbers what it is sent. */
  sent: number[] = [];
  report: DomainBatchReportLike | null = null;
  /** By source row, the key of the row it landed on or repeats. */
  readonly keys = new Map<number, unknown>();
  private _reported = new Map<number, {status: string}>();

  constructor(readonly table: DomainTable, readonly source: DG.DataFrame, readonly mapping: [string, string][],
    count = source.rowCount) {
    this.rows = Array.from({length: count}, () => ({}));
    for (const [target, name] of mapping) {
      const column = source.columns.byName(name);
      for (let i = 0; i < count; i++)
        this.rows[i][target] = column.get(i);
    }
  }

  /** What the caller may write, less what a derivation writes. */
  static targets(table: DomainTable): IProperty[] {
    // an import inserts: the draft policy, so a key typed once on insert is a target too
    const writable = table.access.columnPolicy().forDraft();
    return table.properties.filter((p) =>
      writable.field(p.name!) === 'editable' && !table.derived.targets.includes(p.name!));
  }

  /** What a source column is for, read the way a person reads it: "Chemical name" is the name
   * column and "CAS number" the cas one, so the comparison drops case and everything that is not
   * a letter or a digit, and an exact hit anywhere beats a partial hit somewhere earlier. A column
   * no name matches goes to the one target of its semType. */
  static autoMatch(targets: IProperty[], column: DG.Column): string | undefined {
    const want = DomainTable.matchKey(column.name);
    const names = (p: IProperty) =>
      [p.name ?? '', p.friendlyName ?? ''].map(DomainTable.matchKey).filter((n) => n !== '');
    const hit = (test: (a: string) => boolean) =>
      want === '' ? undefined : targets.find((p) => names(p).some(test))?.name;
    const typed = targets.filter((p) => !isEmpty(column.semType) && p.semType === column.semType);
    return hit((a) => a === want) ?? hit((a) => a.startsWith(want) || want.startsWith(a)) ??
      hit((a) => a.includes(want) || want.includes(a)) ?? (typed.length === 1 ? typed[0].name : undefined);
  }

  /** A column as a form labels it: the declared caption, else the name made presentable. */
  static caption(table: DomainTable, column: string): string {
    const name = table.properties.find((p) => p.name === column)?.friendlyName ?? column.replace(/_/g, ' ');
    return `${name.charAt(0).toUpperCase()}${name.slice(1)}`;
  }

  get keyColumn(): string {
    return DomainPick.nameColumn(this.table.info);
  }

  /** The rows as the table stores them, in place; a row any check failed on is held back. */
  async prepare(): Promise<void> {
    const validators = this.table.validators;
    const problems = [...await this._resolveRefs(), ...await this.table.deriveRows(this.rows),
      ...this.rows.flatMap((row, index) => validators.columns.filter((column) => column in row)
        .map((column) => ({index, column, problem: validators.check(column, row[column], row as RowView)})))
        .filter((p): p is RowProblem => p.problem !== null)];
    const skipped = new Set(problems.filter((p) => p.skip).map((p) => p.index));
    this.problems = problems.filter((p) => p.skip || !skipped.has(p.index)).sort((a, b) => a.index - b.index);
    const held = new Set(this.problems.map((p) => p.index));
    this.sent = this.rows.map((_, i) => i).filter((i) => !held.has(i));
  }

  /** The dry run over the rows to be sent, the held ones added as skips and errors; null where there is none. */
  async validate(options: DomainBatchOptionsLike): Promise<DomainBatchValidationLike | null> {
    const table = this.table.table;
    if (table.validate === undefined)
      return null;
    const report = await table.validate(this.sent.map((i) => this.rows[i]), options);
    const held = this._held();
    const errors = held.filter((row) => row.status === 'error').length;
    return {...report, willSkip: report.willSkip + held.length - errors, errorCount: report.errorCount + errors,
      rows: [...report.rows.map((row) => ({...row, index: this.sent[row.index]})), ...held.map((row) =>
        ({index: row.index, predicted: row.status === 'error' ? 'error' as const : 'skip' as const}))]};
  }

  /** One batch of the rows not held back — none where all or nothing meets an error; reported by source row. */
  async commit(options: DomainBatchOptionsLike): Promise<DomainBatchReportLike> {
    await this.prepare();
    const held = this._held();
    const errors = held.filter((row) => row.status === 'error').length;
    const abort = options.allOrNothing !== false && errors > 0;
    const report: DomainBatchReportLike = abort ? {error: 'validation', rows: []} :
      await this.table.table.batch!(this.sent.map((i) => this.rows[i]), options);
    this.report = {...report, skipped: (report.skipped ?? 0) + held.length - errors,
      errorCount: (report.errorCount ?? report.rows.filter((row) => row.errors?.length).length) + errors,
      rows: [...report.rows.map((row) => ({...row, index: this.sent[row.index]})), ...held]};
    this._reported = new Map(this.report.rows.map((row) => [row.index, row]));
    await this._readKeys();
    return this.report;
  }

  /** A row the report does not name was sent with the rest: the server names at most the first thousand. */
  result(index: number): string {
    const row = this._reported.get(index);
    return row !== undefined ? OUTCOMES[row.status] ?? row.status :
      this.report!.error !== undefined ? 'Not imported' : 'Imported';
  }

  /** The key in the column's display format: DG-000123. */
  keyText(index: number): string {
    const key = this.keys.get(index);
    return isEmpty(key) ? '' : DomainTable.captionText(this.table.properties, this.keyColumn, key);
  }

  /** The rows that landed, as a query by id. */
  shown(): {query: FilterGroup, count: number} | null {
    const ids = this.report!.rows.filter((row) => LANDED.includes(row.status)).map((row) => row.id)
      .filter((id): id is string => !isEmpty(id));
    return ids.length === 0 ? null : {count: ids.length, query: Filters.group('and', [Filters.cond('id', 'in', ids)])};
  }

  /** The source rows beside what became of each: the result, why, and the key. */
  reportFrame(): DG.DataFrame {
    const messages = new Map<number, string[]>();
    const add = (index: number, column: string | undefined, message: string) => messages.set(index,
      [...messages.get(index) ?? [], column ? `${ImportRun.caption(this.table, column)}: ${message}` : message]);
    for (const p of this.problems)
      add(p.index, p.column, p.problem);
    for (const row of this.report!.rows) {
      for (const error of row.errors ?? [])
        add(row.index, error.column, error.message);
    }
    const frame = this.source.clone();
    frame.name = `${this.source.name} import`;
    for (const [name, value] of [[IMPORT_COLUMNS.result, (i: number) => this.result(i)],
      [IMPORT_COLUMNS.message, (i: number) => (messages.get(i) ?? []).join('; ')],
      [ImportRun.caption(this.table, this.keyColumn), (i: number) => this.keyText(i)]] as const) {
      if (frame.columns.contains(name))
        frame.columns.remove(name);
      frame.columns.addNewString(name).init(value);
    }
    return frame;
  }

  /** The rows held back, one line each, as the report lists them. */
  private _held(): {index: number, id: null, status: 'skipped' | 'error'}[] {
    return [...new Map(this.problems.map((p) => [p.index, p.skip ? 'skipped' as const : 'error' as const]))]
      .map(([index, status]) => ({index, id: null, status}));
  }

  /** The key of each row reported with an id, read back: the server may have numbered it. */
  private async _readKeys(): Promise<void> {
    const column = this.keyColumn;
    const idOf = new Map(this.report!.rows.flatMap((row) => row.id === null ? [] : [[row.index, row.id] as const]));
    const ids = [...new Set(idOf.values())];
    const rows = column === 'id' || ids.length === 0 ? [] :
      await this.table.table.query({filter: {property: 'id', operator: '=', value: ids}, columns: [column],
        limit: ids.length});
    const keyOf = new Map(rows.map((row) => [row.id, row[column]]));
    for (const [index, id] of idOf)
      this.keys.set(index, column === 'id' ? id : keyOf.get(id));
  }

  /** A reference column holds the target row's id where a file holds its name ("DGB-000045"): the
   * distinct values of every mapped one looked up by the target's name column, a uuid taken as it
   * stands. A value no row answers is an issue of each row holding it. */
  private async _resolveRefs(): Promise<RowProblem[]> {
    const refs = this.table.properties.filter((p) =>
      DomainPick.isRowRef(p) && this.mapping.some(([target]) => target === p.name));
    const issues: RowProblem[] = [];
    await Promise.all(refs.map(async (prop) => {
      const column = prop.name!;
      const names = new Set(this.rows.map((row) => text(row[column]))
        .filter((v) => v !== '' && !DomainAddress.ID.test(v)));
      const ids = await DomainPick.lookUp(prop.semType!, [...names]);
      for (const [index, row] of this.rows.entries()) {
        const name = text(row[column]);
        if (!names.has(name))
          continue;
        const id = ids.get(name);
        if (id === undefined)
          issues.push({index, column, problem: `No match for "${name}"`});
        else
          row[column] = id;
      }
    }));
    return issues;
  }
}

class ImportFlow {
  private readonly _frame = signal<DG.DataFrame | null>(null);
  /** The source columns, and each one's target choice — kept here, not in the mapping table's cells: the table
   * recycles its rows, and a choice that lived in a cell would be rebuilt, and lose its value, on every scroll. */
  private readonly _columns = signal<string[]>([]);
  private readonly _choices = new Map<string, ChoiceInput>();
  /** Bumped by every edit the gates read — the source, the mode and every choice. */
  private readonly _recheck = signal(0);
  private readonly _report = div([], 'u2-domain-import-report');
  /** Absent where `support.batch` refuses the option. */
  private _mode?: ChoiceInput;
  private _allOrNothing?: BoolInput;
  private _errorOnDuplicate?: BoolInput;
  private _wizard!: Wizard;
  private _committed = false;
  private _result: DomainBatchReportLike | null = null;
  private _settled = false;
  /** The dry run and what it was run for — the mapping, the mode and the flags; the preview posts
   * again only when that key changed. */
  private _checked: {key: string,
    checked: Promise<{report: DomainBatchValidationLike | null, problems: RowProblem[]}>} | undefined;
  /** Bumped by every post: an older dry run answering last must not paint over a newer one. */
  private _previewGen = 0;

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
        {id: 'preview', title: 'Preview', content: () => this._previewStep(), canProceed: () => this._gate()},
        {id: 'report', title: 'Report', content: () => this._report, done: true,
          onActivate: () => void this._commit()},
      ],
      onFinish: () => this._settle(this._result),
      onCancel: () => this._settle(null),
    });
    this._wizard.root.classList.add('u2-domain-import');
    this._wizard.openInDialog(`Import ${DomainSource.noun(this._table.info)}`, {width: 620, height: 520});
  }

  /** The first blocking problem, which is what the NEXT button shows and gates on. */
  private _gate(): string | null {
    return this._problems()[0] ?? null;
  }

  private _settle(report: DomainBatchReportLike | null): void {
    if (this._settled)
      return;
    this._settled = true;
    this._clearChoices();
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
    const form = new Form().addAll([picker, this._mode, this._allOrNothing, this._errorOnDuplicate]
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
    const table = new DataTable<string>({
      columns: [
        {name: 'column', header: 'Source column', render: (column) => column},
        {name: 'target', header: 'Target column', render: (column) => this._choices.get(column)!.root},
      ],
      keyOf: (column) => column,
      labelOf: (column) => column,
      rowHeight: 32,
      items: this._columns,
    });
    table.root.classList.add('u2-domain-import-mapping');
    return table.root;
  }

  /** The server's own verdicts over the mapped cells of the first {@link PREVIEW_ROWS} rows: the
   * dry run runs every check the commit runs, inside a transaction it rolls back, so the preview
   * and the report cannot disagree. Posted once per entry into the step, and again only when the
   * mapping, the mode or the source changed since. */
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
      if (frame !== null && problems.length === 0)
        void this._preview(control, host, frame);
    });
    return host;
  }

  /** The rows as they would land, under the target captions — a structure drawn; the verdict column and the cell
   * marks exist only where there was a dry run to say so. */
  private _rowsTable(owner: Control, frame: DG.DataFrame, mapping: [string, string][],
    verdicts: Map<number, DomainBatchValidationRowLike> | null): HTMLElement {
    const errorOf = (index: number, target: string) => (verdicts?.get(index)?.errors ?? [])
      .find((error) => (error.column ?? '') === (target === VERDICT ? '' : target))?.message ?? null;
    const molecules = new Set(mapping.map(([target]) => target)
      .filter((target) => this._table.properties.find((p) => p.name === target)?.semType === 'Molecule'));
    const table = owner.runInScope(() => new DataTable<number>({
      rowHeight: molecules.size > 0 ? 44 : 24,
      keyOf: (i) => String(i),
      columns: [...(verdicts === null ? [] : [{name: VERDICT, header: 'Result', width: '96px',
        render: (i: number) => {
          const verdict = VERDICTS[verdicts.get(i)?.predicted ?? ''];
          return verdict === undefined ? '' : badge(verdict.label, {variant: verdict.variant});
        }}]),
      ...mapping.map(([target, name]) => {
        const column = frame.columns.byName(name);
        return {name: target, header: ImportRun.caption(this._table, target), render: (i: number) => {
          const value = column.get(i);
          return molecules.has(target) && !isEmpty(value) ? THUMBNAIL.listItem!(String(value)) : text(value);
        }};
      })],
      cellState: {
        isChanged: () => false,
        errorOf: (key, target) => {
          const message = errorOf(Number(key), target);
          return message === null ? null : {message, kind: 'error'};
        },
      },
    }));
    table.setItems(Array.from({length: Math.min(frame.rowCount, PREVIEW_SHOWN)}, (_, i) => i));
    table.root.classList.add('u2-domain-import-rows');
    return table.root;
  }

  /** The rows prepared as the commit prepares them, and the dry run over them where the storage
   * has one; without it, the mapped rows alone. */
  private async _preview(owner: Control, host: HTMLElement, frame: DG.DataFrame): Promise<void> {
    const mapping = [...this._mapped()];
    const options = this._batchOptions();
    const scanned = Math.min(frame.rowCount, PREVIEW_ROWS);
    const key = JSON.stringify([mapping, options]);
    const gen = ++this._previewGen;
    let report: DomainBatchValidationLike | null;
    let problems: RowProblem[];
    try {
      if (this._checked?.key !== key) {
        const run = new ImportRun(this._table, frame, mapping, scanned);
        this._checked = {key, checked: (async () => {
          await run.prepare();
          return {problems: run.problems, report: await run.validate(options)};
        })()};
        host.replaceChildren(span(`Checking ${scanned} rows…`, 'u2-domain-import-summary'));
      }
      ({report, problems} = await this._checked.checked);
      if (gen !== this._previewGen)
        return;
    } catch (e) {
      if (gen !== this._previewGen)
        return;
      // the commit is the authority: a preview that cannot reach the server says so and lets NEXT
      // through
      this._checked = undefined;
      host.replaceChildren(span(DomainErrors.message(e), 'u2-domain-import-problem'));
      return;
    }
    if (owner.scope.isDisposed)
      return;
    host.replaceChildren(this._rowsTable(owner, frame, mapping,
      report === null ? null : new Map(report.rows.map((row) => [row.index, row]))));
    if (report === null) {
      const shown = Math.min(frame.rowCount, PREVIEW_SHOWN);
      const lead = shown < frame.rowCount ? `Showing the first ${shown.toLocaleString()} of ` +
        `${frame.rowCount.toLocaleString()} rows — rows` : 'Rows';
      host.append(span(`${lead} are checked when imported.`, 'u2-domain-import-summary'));
    }
    else {
      host.append(span(`${report.willInsert} will be added, ${report.willUpdate} updated, ` +
        `${report.willSkip} skipped, ${plural(report.errorCount, 'row has errors', 'rows have errors')}.`,
      'u2-domain-import-summary'));
      // the counts are what the batch WOULD do row by row; all-or-nothing makes one bad row the
      // verdict on all of them, and the report would otherwise contradict this line
      if (options.allOrNothing !== false && report.errorCount > 0) {
        host.append(span(this._allOrNothing === undefined ?
          'This import is all or nothing — nothing will be imported.' :
          'Nothing will be imported while "All or nothing" is on.', 'u2-domain-import-problem'));
      }
      if (scanned < frame.rowCount) {
        host.append(span(`Checked the first ${scanned.toLocaleString()} of ` +
          `${frame.rowCount.toLocaleString()} rows — the rest are checked on import.`,
        'u2-domain-import-summary'));
      }
    }
    host.append(...this._issues(owner, ImportFlow._collect(problems, report?.rows ?? [])));
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

  /** One target choice per source column, set to what it auto-matches; the asterisk the forms use in the option
   * itself — a picker's own list is all the user sees — and on the choice standing on a required target. */
  private _remap(frame: DG.DataFrame | null): void {
    const required = (target: string | null) => this._targets.some((p) => p.name === target && p.nullable === false);
    const items: ChoiceItem[] = this._targets.map((p) =>
      ({value: p.name!, label: `${ImportRun.caption(this._table, p.name!)}${p.nullable === false ? ' *' : ''}`}));
    this._clearChoices();
    for (const column of frame?.columns.toList() ?? []) {
      // owned by hand, not by a scope: a remap may run before the wizard that shows it is built
      const choice = new ChoiceInput({items: [SKIP, ...items], inline: true, name: column.name, nullable: false,
        value: ImportRun.autoMatch(this._targets, column) ?? SKIP});
      choice.effect(() => choice.root.classList.toggle('u2-input-required', required(choice.value.value)));
      choice.value.subscribe(() => this._touch());
      this._choices.set(column.name, choice);
    }
    this._columns.value = frame?.columns.names() ?? [];
    this._touch();
  }

  private _clearChoices(): void {
    for (const choice of this._choices.values())
      choice.dispose();
    this._choices.clear();
  }

  /** Target → source column, in mapping order; a target claimed twice keeps the first. */
  private _mapped(): Map<string, string> {
    const out = new Map<string, string>();
    for (const [column, choice] of this._choices) {
      const target = choice.value.peek();
      if (target !== null && target !== SKIP && !out.has(target))
        out.set(target, column);
    }
    return out;
  }

  /** The targets more than one source column is mapped to. */
  private _duplicates(): string[] {
    const claimed = new Set<string>();
    const twice = new Set<string>();
    for (const choice of this._choices.values()) {
      const target = choice.value.peek();
      if (target === null || target === SKIP)
        continue;
      if (claimed.has(target))
        twice.add(target);
      claimed.add(target);
    }
    return [...twice];
  }

  /** Everything that blocks the import, in the order `domain_import_dialog.dart` reports it. */
  private _problems(): string[] {
    this._recheck.value;
    if (this._frame.value === null)
      return ['Choose an open table or a file'];
    const problems = this._duplicates().map((target) =>
      `Two source columns are mapped to "${target}" — keep one and skip the other.`);
    const mapped = this._mapped();
    if (mapped.size === 0)
      problems.push('Map at least one column.');
    if (this._mode?.value.peek() === 'upsert') {
      for (const key of this._table.info.businessKey) {
        if (!mapped.has(key)) {
          problems.push('Upsert merges by the business key — map a column to ' +
            `"${ImportRun.caption(this._table, key)}".`);
        }
      }
    } else {
      for (const p of this._targets) {
        if (p.nullable === false && (p.defaultValue === undefined || p.defaultValue === null) && !mapped.has(p.name!)) {
          problems.push(`Required column "${ImportRun.caption(this._table, p.name!)}" is not mapped — ` +
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
    this._report.replaceChildren(span(`Importing ${frame.rowCount} rows…`, 'u2-domain-import-summary'));
    // reading the frame is part of the import: a column the picked frame no longer has belongs in
    // the report with everything else the server refuses
    let run: ImportRun;
    try {
      run = new ImportRun(this._table, frame, [...this._mapped()]);
      this._result = await run.commit(this._batchOptions());
    } catch (e) {
      this._report.replaceChildren(span(DomainErrors.message(e), 'u2-domain-import-problem'));
      notify.error(DomainErrors.message(e));
      return;
    }
    this._render(run);
  }

  private _render(run: ImportRun): void {
    const report = run.report!;
    const failed = report.error !== undefined;
    // a warehouse refusal carries the failing rows and no totals
    const errorCount = report.errorCount ?? report.rows.filter((row) => row.errors?.length).length;
    this._report.replaceChildren(span(failed ?
      `Import aborted — ${plural(errorCount, 'row has errors', 'rows have errors')}; ` +
      'nothing was committed.' :
      `${ImportFlow._landed(report)}, ${report.skipped ?? 0} skipped, ${errorCount} failed.`,
    failed ? 'u2-domain-import-problem' : 'u2-domain-import-summary'));
    this._report.append(...this._issues(this._wizard, ImportFlow._collect(run.problems, report.rows)),
      ...this._keys(run), this._links(run));
    if (failed) {
      this._report.append(...this._wayBack(errorCount + (report.skipped ?? 0)));
      notify.error('Import failed — nothing was committed. See the report.');
    }
    else
      notify.info(`Imported: ${ImportFlow._landed(report)}`);
  }

  /** The rows that landed or repeat one, by the key the rest of the app names them with. */
  private _keys(run: ImportRun): HTMLElement[] {
    const rows = [...run.keys.keys()].sort((a, b) => a - b);
    if (rows.length === 0)
      return [];
    const table = this._wizard.runInScope(() => new DataTable<number>({rowHeight: 24, keyOf: (i) => String(i),
      columns: [
        {name: 'row', header: 'Source row', width: '90px', render: (i) => String(i + 1)},
        {name: 'result', header: 'Result', width: '96px', render: (i) => run.result(i)},
        {name: 'key', header: ImportRun.caption(this._table, run.keyColumn), render: (i) => run.keyText(i)},
      ]}));
    table.setItems(rows);
    table.root.classList.add('u2-domain-import-rows');
    return [table.root];
  }

  /** "Show the N imported rows", where the caller can show them, and every row's outcome as a table. */
  private _links(run: ImportRun): HTMLElement {
    const shown = run.shown();
    const show = this._options.show;
    return div([...(shown === null || show === undefined ? [] :
      [link(`Show the ${plural(shown.count, 'imported row', 'imported rows')}`, () => {
        this._wizard.next();
        show(shown.query);
      })]), link('Open as table', () => grok.shell.addTableView(run.reportFrame()))], 'u2-domain-import-links');
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
  private _wayBack(notImported: number): HTMLElement[] {
    const valid = (this._frame.peek()?.rowCount ?? 0) - notImported;
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

  /** The lines a preview or a report lists — the rows' own failures, then the backend's row by
   * row — at most {@link ISSUE_CAP} of them, and how many there are in all. */
  private static _collect(problems: RowProblem[], rows: {index: number, status?: string,
    errors?: {column?: string, message: string}[]}[]): {issues: RowProblem[], total: number} {
    const issues = problems.slice(0, ISSUE_CAP);
    let total = problems.length;
    const add = (issue: RowProblem) => {
      total++;
      if (issues.length < ISSUE_CAP)
        issues.push(issue);
    };
    // a cell the import explained already (a name no row holds) needs no second word from the backend
    const explained = new Set(problems.map((p) => `${p.index}|${p.column}`));
    for (const row of rows) {
      const errors = row.errors ?? [];
      for (const error of errors.filter((e) => !explained.has(`${row.index}|${e.column ?? ''}`)))
        add({index: row.index, column: error.column ?? '', problem: error.message});
      if (errors.length === 0 && row.status === 'duplicate')
        add({index: row.index, column: '', problem: 'Duplicate business key — skipped'});
    }
    return {issues, total};
  }

  private _issues(owner: Control, {issues, total}: {issues: RowProblem[], total: number}): HTMLElement[] {
    if (issues.length === 0)
      return [];
    // the wire counts rows from zero; the person reading this counts the source frame from one
    const table = owner.runInScope(() => new BasicTable<RowProblem>({items: issues, columns: [
      {header: 'Source row', render: (x) => String(x.index + 1), width: '90px'},
      {header: 'Column', render: (x) => ImportRun.caption(this._table, x.column), width: '140px'},
      {header: 'Message', render: (x) => x.problem},
    ]}));
    const host = div([table.root], 'u2-domain-import-issues');
    if (total > issues.length)
      host.append(span(`…and ${total - issues.length} more`));
    return [host];
  }

}

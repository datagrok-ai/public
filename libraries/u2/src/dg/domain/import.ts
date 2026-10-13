/* `domains.import` — the import wizard (WO 3-9), the u2 half of `domain_import_dialog.dart`:
   **source** (any open frame, or any file the platform can read, through `ui.input.table`),
   **mapping** (a `DataTable` of source column → target, auto-matched by name; only the columns
   the caller may write are offered), **preview** (the server's dry run over the first 1000 rows —
   `validate`, which is the commit's own checks inside a rolled-back transaction; the mapped rows
   alone where the storage has no dry run) and **report** (the commit's). What is posted is the
   MAPPED columns only, under their target names, so a skipped or renamed source column never
   reaches the server and the source frame is never touched. The options offered are the ones
   `support.batch` declares: an option the storage refuses is not on the form and not sent. */
import * as ui from 'datagrok-api/ui';
import type * as DG from 'datagrok-api/dg';
import {signal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {button, div, divV, span} from '../../core/elements.js';
import {plural, text} from '../../core/text.js';
import type {IProperty} from '../../core/property-like.js';
import {Wizard} from '../../components/containers/wizard.js';
import {DataTable} from '../../components/collections/data-table.js';
import {BasicTable} from '../../components/collections/table.js';
import {ChoiceInput} from '../../components/inputs/choice-input.js';
import type {ChoiceItem} from '../../components/inputs/choice-input.js';
import {BoolInput} from '../../components/inputs/bool-input.js';
import {badge} from '../../components/display/badge.js';
import type {BadgeVariant} from '../../components/display/badge.js';
import {Form} from '../../components/forms/form.js';
import {notify} from '../../components/display/notify.js';
import type {DomainBatchOptionsLike, DomainBatchReportLike, DomainBatchValidationLike,
  DomainBatchValidationRowLike} from '../../sources/domain-backend.js';
import {fromDartInput} from '../inputs/from-dart-input.js';
import {DomainTable} from './index.js';
import {DomainErrors} from './errors.js';

export interface DomainImportOptions {
  /** The frame to import, when the caller already has one; the picker still offers the rest. */
  source?: DG.DataFrame;
}

/** What a source column maps to when it maps to nothing. */
const SKIP = '(skip)';
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
  const targets = table.properties.filter((p) => writable.field(p.name!) === 'editable');
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
  private _wizard!: Wizard;
  private _committed = false;
  private _result: DomainBatchReportLike | null = null;
  private _settled = false;
  /** The dry run and what it was run for — the mapping, the mode and the flags; the preview posts
   * again only when that key changed. */
  private _checked: {key: string, report: Promise<DomainBatchValidationLike>} | undefined;
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
    this._wizard.openInDialog(`Import ${ImportFlow._plural(this._table)}`, {width: 620, height: 520});
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
      if (frame === null || problems.length > 0)
        return;
      if (this._table.table.validate === undefined) {
        const shown = Math.min(frame.rowCount, PREVIEW_SHOWN);
        const lead = shown < frame.rowCount ? `Showing the first ${shown.toLocaleString()} of ` +
          `${frame.rowCount.toLocaleString()} rows — rows` : 'Rows';
        host.append(this._rowsTable(control, frame, [...this._mapping()], null),
          span(`${lead} are checked when imported.`, 'u2-domain-import-summary'));
        return;
      }
      void this._preview(control, host, frame);
    });
    return host;
  }

  /** The rows as they would land, under the target captions; the verdict column and the cell
   * marks exist only where there was a dry run to say so. */
  private _rowsTable(owner: Control, frame: DG.DataFrame, mapping: [string, string][],
    verdicts: Map<number, DomainBatchValidationRowLike> | null): HTMLElement {
    const errorOf = (index: number, target: string) => (verdicts?.get(index)?.errors ?? [])
      .find((error) => (error.column ?? '') === (target === VERDICT ? '' : target))?.message ?? null;
    const table = owner.runInScope(() => new DataTable<number>({
      rowHeight: 24,
      keyOf: (i) => String(i),
      columns: [...(verdicts === null ? [] : [{name: VERDICT, header: 'Result', width: '96px',
        render: (i: number) => {
          const verdict = VERDICTS[verdicts.get(i)?.predicted ?? ''];
          return verdict === undefined ? '' : badge(verdict.label, {variant: verdict.variant});
        }}]),
      ...mapping.map(([target, column]) => ({
        name: target, header: ImportFlow._caption(this, target),
        render: (i: number) => text(frame.get(column, i)),
      }))],
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

  private async _preview(owner: Control, host: HTMLElement, frame: DG.DataFrame): Promise<void> {
    const mapping = [...this._mapping()];
    const options = this._batchOptions();
    const scanned = Math.min(frame.rowCount, PREVIEW_ROWS);
    const key = JSON.stringify([mapping, options]);
    const gen = ++this._previewGen;
    let report: DomainBatchValidationLike;
    try {
      if (this._checked?.key !== key) {
        this._checked = {key, report: this._table.table.validate!(this._payload(frame, mapping, scanned), options)};
        host.replaceChildren(span(`Checking ${scanned} rows…`, 'u2-domain-import-summary'));
      }
      report = await this._checked.report;
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
    host.replaceChildren(this._rowsTable(owner, frame, mapping, new Map(report.rows.map((row) => [row.index, row]))));
    host.append(span(`${report.willInsert} will be added, ${report.willUpdate} updated, ` +
      `${report.willSkip} skipped, ${plural(report.errorCount, 'row has errors', 'rows have errors')}.`,
    'u2-domain-import-summary'));
    // the counts are what the batch WOULD do row by row; all-or-nothing makes one bad row the
    // verdict on all of them, and the report would otherwise contradict this line
    if (options.allOrNothing !== false && report.errorCount > 0)
      host.append(span(this._allOrNothing === undefined ? 'This import is all or nothing — nothing will be imported.' :
        'Nothing will be imported while "All or nothing" is on.', 'u2-domain-import-problem'));
    if (scanned < frame.rowCount) {
      host.append(span(`Checked the first ${scanned.toLocaleString()} of ` +
        `${frame.rowCount.toLocaleString()} rows — the rest are checked on import.`,
      'u2-domain-import-summary'));
    }
    const issues: Issue[] = [];
    let total = 0;
    for (const row of report.rows) {
      for (const error of row.errors ?? []) {
        total++;
        if (issues.length < ISSUE_CAP)
          issues.push({row: row.index, column: error.column ?? '', message: error.message});
      }
    }
    if (issues.length > 0)
      host.append(ImportFlow._issues(owner, issues, total));
  }

  /** The mapped rows, under their target names — what is posted, and the only thing that is. */
  private _payload(frame: DG.DataFrame, mapping: [string, string][], count: number):
    Record<string, unknown>[] {
    const rows: Record<string, unknown>[] = [];
    for (let i = 0; i < count; i++) {
      const row: Record<string, unknown> = {};
      for (const [target, column] of mapping)
        row[target] = frame.get(column, i);
      rows.push(row);
    }
    return rows;
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

  /** One target choice per source column, auto-matched case-insensitively against the property's
   * name and its caption. The choices are kept here, not in the mapping table's cells: the table
   * recycles its rows, and a choice that lived in a cell would be rebuilt — and lose its value —
   * on every scroll. */
  private _remap(frame: DG.DataFrame | null): void {
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
      const choice = new ChoiceInput({items, value: this._autoMatch(column), inline: true, name: column,
        nullable: false});
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
    const prop = flow._targets.find((p) => p.name === column);
    return prop === undefined ? column : ImportFlow.captionOf(prop);
  }

  /** The table's plural name, as the user reads it. */
  private static _plural(table: DomainTable): string {
    return (table.info.pluralName || 'rows').replace(/_/g, ' ').toLowerCase();
  }

  /** What a source column is for, read the way a person reads it: "Chemical name" is the name
   * column and "CAS number" the cas one, so the comparison drops case and everything that is not
   * a letter or a digit, and an exact hit anywhere beats a partial hit somewhere earlier. */
  private _autoMatch(column: string): string {
    const plain = (s: string) => s.toLowerCase().replace(/[^a-z0-9]/g, '');
    const want = plain(column);
    if (want === '')
      return SKIP;
    const names = (p: IProperty) => [plain(p.name ?? ''), plain(p.friendlyName ?? '')].filter((n) => n !== '');
    const hit = (test: (a: string) => boolean) => this._targets.find((p) => names(p).some(test))?.name;
    return hit((a) => a === want) ?? hit((a) => a.startsWith(want) || want.startsWith(a)) ??
      hit((a) => a.includes(want) || want.includes(a)) ?? SKIP;
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
    this._report.replaceChildren(span(`Importing ${frame.rowCount} rows…`, 'u2-domain-import-summary'));
    // reading the frame is part of the import: a column the picked frame no longer has belongs in
    // the report with everything else the server refuses
    try {
      this._result = await this._table.table.batch!(
        this._payload(frame, [...this._mapping()], frame.rowCount), this._batchOptions());
    } catch (e) {
      this._report.replaceChildren(span(DomainErrors.message(e), 'u2-domain-import-problem'));
      notify.error(DomainErrors.message(e));
      return;
    }
    this._render(this._result);
  }

  private _render(report: DomainBatchReportLike): void {
    const failed = report.error !== undefined;
    // a warehouse refusal carries the failing rows and no totals
    const errorCount = report.errorCount ?? report.rows.filter((row) => row.errors?.length).length;
    this._report.replaceChildren(span(failed ?
      `Import aborted — ${plural(errorCount, 'row has errors', 'rows have errors')}; ` +
      'nothing was committed.' :
      `${ImportFlow._landed(report)}, ${report.skipped ?? 0} skipped, ${errorCount} failed.`,
    failed ? 'u2-domain-import-problem' : 'u2-domain-import-summary'));
    const issues: Issue[] = [];
    let total = 0;
    for (const row of report.rows) {
      const errors = row.errors ?? [];
      for (const error of errors) {
        total++;
        if (issues.length < ISSUE_CAP)
          issues.push({row: row.index, column: error.column ?? '', message: error.message});
      }
      if (errors.length === 0 && row.status === 'duplicate') {
        total++;
        if (issues.length < ISSUE_CAP)
          issues.push({row: row.index, column: '', message: 'Duplicate business key — skipped'});
      }
    }
    if (issues.length > 0)
      this._report.append(ImportFlow._issues(this._wizard, issues, total));
    if (failed) {
      this._report.append(...this._wayBack(errorCount));
      notify.error('Import failed — nothing was committed. See the report.');
    }
    else
      notify.info(`Imported: ${ImportFlow._landed(report)}`);
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
  private _wayBack(errorCount: number): HTMLElement[] {
    const valid = (this._frame.peek()?.rowCount ?? 0) - errorCount;
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

  private static _issues(owner: Control, issues: Issue[], total: number): HTMLElement {
    // the wire counts rows from zero; the person reading this counts the source frame from one
    const table = owner.runInScope(() => new BasicTable<Issue>({items: issues, columns: [
      {header: 'Source row', render: (x) => String(x.row + 1), width: '90px'},
      {header: 'Column', render: (x) => x.column, width: '140px'},
      {header: 'Message', render: (x) => x.message},
    ]}));
    const host = div([table.root], 'u2-domain-import-issues');
    if (total > issues.length)
      host.append(span(`…and ${total - issues.length} more`));
    return host;
  }

}

interface Issue {row: number; column: string; message: string}

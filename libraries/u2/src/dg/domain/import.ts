/* `domains.import` — the import wizard (WO 3-9), the u2 half of `domain_import_dialog.dart`:
   **source** (any open frame, or any file the platform can read, through `ui.input.table`),
   **mapping** (a `DataTable` of source column → target, auto-matched by name; only the columns
   the caller may write are offered), **preview** (the schema's own rules over the first 10k rows
   — advisory, since the batch endpoint has no dry run) and **report** (the server's, which is the
   authority). What is posted is the MAPPED columns only, under their target names, so a skipped
   or renamed source column never reaches the server and the source frame is never touched. */
import * as ui from 'datagrok-api/ui';
import type * as DG from 'datagrok-api/dg';
import {signal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {Filters} from '../../core/filter/index.js';
import {div, divV, span} from '../../core/elements.js';
import {text} from '../../core/text.js';
import type {IProperty} from '../../core/property-like.js';
import {Wizard} from '../../components/containers/wizard.js';
import {DataTable} from '../../components/collections/data-table.js';
import {BasicTable} from '../../components/collections/table.js';
import {ChoiceInput} from '../../components/inputs/choice-input.js';
import type {ChoiceItem} from '../../components/inputs/choice-input.js';
import {BoolInput} from '../../components/inputs/bool-input.js';
import {Form} from '../../components/forms/form.js';
import {notify} from '../../components/display/notify.js';
import {MemoryEditState} from '../../sources/edit-state.js';
import type {DomainBatchOptionsLike, DomainBatchReportLike,
  DomainTransactionOpLike} from '../../sources/domain-backend.js';
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
/** How many business keys one upsert lookup asks about at a time (a backend with no `/batch`). */
const KEY_CHUNK = 100;
/** Issue lines rendered in a preview or a report; the counts stay exact. */
const ISSUE_CAP = 20;
/** Rows the preview draws — what the import looks like, not all of it. */
const PREVIEW_SHOWN = 20;

/** One row of the mapping step — the target is the choice input the flow keeps for that column. */
interface Mapping {column: string}

/** Opens the wizard over `table`; resolves to the server's report, or to null where it was
 * cancelled. Exported as `domains.import` — `import` cannot name a function. */
export function openImport(table: DomainTable, options: DomainImportOptions = {}):
  Promise<DomainBatchReportLike | null> {
  // column security alone: on a row-mode table the table-level `insert`/`edit` are false negatives
  // (GOAL "Access"), and the server refuses a column the caller may not write anyway
  const writable = table.access.narrow({edit: true, insert: true});
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
  private _mode!: ChoiceInput;
  private _allOrNothing!: BoolInput;
  private _errorOnDuplicate!: BoolInput;
  private _wizard!: Wizard;
  private _committed = false;
  private _result: DomainBatchReportLike | null = null;
  private _settled = false;

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
    this._mode = new ChoiceInput({label: 'Mode', name: 'mode', value: 'insert', nullable: false,
      items: [{value: 'insert', label: 'Add new rows'}, {value: 'upsert', label: 'Add or update'}]});
    this._allOrNothing = new BoolInput({label: 'All or nothing', value: true, name: 'allOrNothing',
      tooltipText: 'Any invalid row aborts the whole import'});
    this._errorOnDuplicate = new BoolInput({label: 'Error on duplicate', value: false, name: 'errorOnDuplicate',
      tooltipText: 'Treat business-key duplicates as errors instead of skips (insert mode)'});
    const form = new Form().addAll([picker, this._mode, this._allOrNothing, this._errorOnDuplicate]);
    form.effect(() => {
      const frame = picker.value.value;
      summary.textContent = frame === null ? '' :
        `${frame.name}: ${frame.rowCount} rows, ${frame.columns.names().length} columns`;
      this._frame.value = frame;
      this._remap(frame);
      this._touch();
    });
    form.effect(() => {
      this._errorOnDuplicate.enabled = this._mode.value.value === 'insert';
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

  /** The client-side pass: the schema's own rules over the mapped cells of the first
   * {@link PREVIEW_ROWS} rows. Advisory — the server re-validates everything it is sent. */
  private _previewStep(): HTMLElement {
    const host = div([], 'u2-domain-import-preview');
    const control = new Control(host);
    control.effect(() => {
      this._recheck.value;
      const frame = this._frame.value;
      const problems = this._problems();
      host.replaceChildren(...problems.map((p) => span(p, 'u2-domain-import-problem')));
      if (frame === null || problems.length > 0)
        return;
      const scanned = Math.min(frame.rowCount, PREVIEW_ROWS);
      const issues: Issue[] = [];
      const bad = new Set<number>();
      let total = 0;
      for (const [target, column] of this._mapping()) {
        const prop = this._targets.find((p) => p.name === target)!;
        for (let i = 0; i < scanned; i++) {
          const message = ImportFlow.problemOf(prop, frame.get(column, i));
          if (message === null)
            continue;
          total++;
          bad.add(i);
          if (issues.length < ISSUE_CAP)
            issues.push({row: i, column: target, message});
        }
      }
      // the rows as they would land, under the target captions, with every refused cell marked:
      // a preview that shows nothing is only a promise that something was checked
      const mapping = [...this._mapping()];
      const shown = Array.from({length: Math.min(frame.rowCount, PREVIEW_SHOWN)}, (_, i) => i);
      const table = control.runInScope(() => new DataTable<number>({
        rowHeight: 24,
        keyOf: (i) => String(i),
        columns: mapping.map(([target, column]) => ({
          name: target, header: ImportFlow._caption(this, target),
          render: (i: number) => text(frame.get(column, i)),
        })),
        cellState: {
          isChanged: () => false,
          errorOf: (key, target) => {
            const [, column] = mapping.find(([t]) => t === target)!;
            const prop = this._targets.find((p) => p.name === target)!;
            const message = ImportFlow.problemOf(prop, frame.get(column, Number(key)));
            return message === null ? null : {message, kind: 'error'};
          },
        },
      }));
      table.setItems(shown);
      table.root.classList.add('u2-domain-import-rows');
      host.append(table.root);
      const tail = scanned < frame.rowCount ? ` (first ${scanned} rows checked)` : '';
      host.append(span(bad.size === 0 ?
        `${frame.rowCount} rows look valid${tail} — the server re-validates on import.` :
        `${bad.size} of ${scanned} checked rows have problems${tail}.`, 'u2-domain-import-summary'));
      if (issues.length > 0)
        host.append(ImportFlow._issues(control, issues, total));
    });
    return host;
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
    if (this._mode.value.peek() === 'upsert') {
      const businessKey = this._table.info.businessKey;
      if (businessKey.length === 0)
        problems.push('Upsert merges by the business key, and this table declares none — import as insert.');
      for (const key of businessKey) {
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
    const mode = this._mode.value.peek() === 'upsert' ? 'upsert' : 'insert';
    // reading the frame is part of the import: a column the picked frame no longer has belongs in
    // the report with everything else the server refuses
    try {
      const mapping = this._mapping();
      const rows: Record<string, unknown>[] = [];
      for (let i = 0; i < frame.rowCount; i++) {
        const row: Record<string, unknown> = {};
        for (const [target, column] of mapping)
          row[target] = frame.get(column, i);
        rows.push(row);
      }
      this._result = await ImportFlow.post(this._table, rows, {mode, allOrNothing: this._allOrNothing.value.peek(),
        errorOnDuplicate: mode === 'insert' && this._errorOnDuplicate.value.peek()});
    } catch (e) {
      this._report.replaceChildren(span(DomainErrors.message(e), 'u2-domain-import-problem'));
      notify.error(DomainErrors.message(e));
      return;
    }
    this._render(this._result);
  }

  private _render(report: DomainBatchReportLike): void {
    const failed = report.error !== undefined;
    this._report.replaceChildren(span(failed ?
      `Import aborted — ${report.errorCount} row(s) with errors; nothing was committed.` :
      `${report.inserted} inserted, ${report.updated} updated, ${report.skipped} skipped, ` +
      `${report.errorCount} failed.`, failed ? 'u2-domain-import-problem' : 'u2-domain-import-summary'));
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
    if (failed)
      notify.error('Import failed — nothing was committed. See the report.');
    else
      notify.info(`Imported: ${report.inserted} inserted, ${report.updated} updated`);
  }

  private static _issues(owner: Control, issues: Issue[], total: number): HTMLElement {
    const table = owner.runInScope(() => new BasicTable<Issue>({items: issues, columns: [
      {header: 'Row', render: (x) => String(x.row), width: '60px'},
      {header: 'Column', render: (x) => x.column, width: '140px'},
      {header: 'Message', render: (x) => x.message},
    ]}));
    const host = div([table.root], 'u2-domain-import-issues');
    if (total > issues.length)
      host.append(span(`…and ${total - issues.length} more`));
    return host;
  }

  /** The rows posted as ONE write: the backend's bulk endpoint where it declares one (the
   * platform's `/batch`, which owns the upsert merge and the per-row report), else one
   * transaction — the memory backend's path, which has no bulk endpoint. */
  static async post(table: DomainTable, rows: Record<string, unknown>[],
    options: DomainBatchOptionsLike): Promise<DomainBatchReportLike> {
    const handle = table.table;
    if (handle.batch !== undefined)
      return handle.batch(rows, options);
    const existing = options.mode === 'upsert' ? await ImportFlow._byKey(table, rows) : new Map<string, string>();
    const ops: DomainTransactionOpLike[] = rows.map((values) => {
      const id = existing.get(ImportFlow._key(table, values));
      return id === undefined ? {op: 'insert' as const, table: table.address, values} :
        {op: 'update' as const, table: table.address, id, values};
    });
    const results = await handle.transaction(ops);
    return {
      inserted: ops.filter((op) => op.op === 'insert').length,
      updated: ops.filter((op) => op.op === 'update').length,
      skipped: 0, errorCount: 0,
      rows: results.map((r, index) => ({index, id: r.id ?? null,
        status: ops[index].op === 'insert' ? 'inserted' : 'updated'})),
    };
  }

  /** Business key → id, for the rows an upsert would merge into — asked for by key, in chunks, so
   * a table with more rows than one page does not silently insert duplicates. */
  private static async _byKey(table: DomainTable, rows: Record<string, unknown>[]): Promise<Map<string, string>> {
    const key = table.info.businessKey;
    const wanted = new Set(rows.map((values) => ImportFlow._key(table, values)));
    const heads = [...new Set(rows.map((values) => String(values[key[0]] ?? '')))];
    const out = new Map<string, string>();
    for (let at = 0; at < heads.length; at += KEY_CHUNK) {
      const filter = Filters.toDomainTree(Filters.group('and',
        [Filters.cond(key[0], 'in', heads.slice(at, at + KEY_CHUNK))]));
      for (const row of await table.table.query({filter, columns: ['id', ...key], limit: KEY_CHUNK * 10})) {
        const spelled = ImportFlow._key(table, row);
        if (wanted.has(spelled))
          out.set(spelled, String(row.id));
      }
    }
    return out;
  }

  private static _key(table: DomainTable, values: Record<string, unknown>): string {
    return table.info.businessKey.map((column) => String(values[column] ?? '')).join('\u0000');
  }

  /** The schema's rules on one cell, plus the coercion the batch engine applies to a TEXT cell:
   * an int column parses the whole string (`'3.0'` is refused, as Dart's `int.parse` refuses it),
   * a bool takes only `true`/`false`, and a datetime is left to the server, whose parse is the
   * authority (`domain_import_dialog.dart` `_previewValidate`). */
  static problemOf(prop: IProperty, value: unknown): string | null {
    const type = prop.propertyType ?? prop.type;
    if (typeof value !== 'string')
      return MemoryEditState.problemOf(prop, value);
    const text = value.trim();
    if (text === '')
      return MemoryEditState.problemOf(prop, null);
    if (type === 'int') {
      return /^[+-]?\d+$/.test(text) ? MemoryEditState.problemOf(prop, Number(text)) :
        `Integer value expected, passed: "${value}"`;
    }
    if (type === 'double' || type === 'float' || type === 'num') {
      const n = Number(text);
      return Number.isNaN(n) ? `Numerical value expected, passed: "${value}"` :
        MemoryEditState.problemOf(prop, n);
    }
    if (type === 'bool') {
      return ['true', 'false'].includes(text.toLowerCase()) ? null :
        `Boolean value expected (true/false), passed: "${value}"`;
    }
    if (type === 'datetime')
      return null;
    return MemoryEditState.problemOf(prop, value);
  }
}

interface Issue {row: number; column: string; message: string}

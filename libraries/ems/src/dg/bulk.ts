/* `domains.bulkEdit` — one value written into many rows at once (WO 3-8), through the platform's
   `POST …/{table}/update`. One aligned form: what the edit applies to, the fields it changes (a
   picker over the columns the caller may write), and an editor for each chosen field only — so
   "clear this field" (chosen, left empty) and "leave it alone" (not chosen) are different requests.
   Under each chosen field a muted line says what it does to the target: the rows' shared value or
   how many differ, or that it will be cleared. The target is the selection, or everything the
   list's filter matches; the server narrows both by the Edit predicate and answers how many rows
   it touched.

   Only targets that can actually be written are offered. Whatever blocks the button is its
   tooltip, and a refusal — the table's validators over the rows as they would land, or the
   server's — is the dialog's error line, never a balloon: a balloon says what landed, nothing else.

   A bulk edit is NOT part of the unit of work — it lands the moment the button is pressed — so the
   session's pending changes are settled through `confirmDiscard` first. */
import {signal, computed, Filters, span, Dialog, Input, MultiSelect, RadioInput, Rows, notify, plural, text}
  from '@datagrok-libraries/u2';
import type {IProperty, ReadonlySignal, RowView} from '@datagrok-libraries/u2';
import {confirmDiscard} from '../core/session.js';
import type {DomainSource} from '../core/domain-source.js';
import type {DomainQueryLike} from '../core/domain-backend.js';
import {ObjectForm, propertyForm} from '@datagrok-libraries/u2/src/dg/index.js';
import {DomainTable, firstProblem} from './facade.js';
import {DomainErrors} from './errors.js';

export interface DomainBulkEditOptions {
  /** The columns offered, in this order; every column the caller may write by default. */
  columns?: string[];
}

/** The server's `maxUpdateWhereRows` — a selection past it is refused by name rather than
 * silently truncated to the first 1000 rows the server happens to order first. */
const MAX_ROWS = 1000;
const SELECTED = 'selected';
const MATCHING = 'matching';

interface BulkTarget {
  kind: typeof SELECTED | typeof MATCHING;
  label: string;
  count: number;
}

/** Writes the chosen columns into the rows the dialog's target names; resolves to how many rows
 * the server touched, or null where the dialog was cancelled. */
export async function bulkEdit(source: DomainSource, options: DomainBulkEditOptions = {}): Promise<number | null> {
  const table = DomainTable.of(source)?.table;
  if (table?.updateWhere === undefined) {
    notify.error(`${source.table}: the backend does not support bulk edits`);
    return null;
  }
  // the trash mode of the source, not a bound on its access: nothing in it is written until the
  // rows are restored, and the column policy below reads column security alone
  if (source.readOnly.peek()) {
    notify.warning(`${source.table}: deleted rows are read-only until they are restored`);
    return null;
  }
  const writable = source.access.peek().columnPolicy();
  const properties = source.schema.properties;
  const props: IProperty[] = [];
  for (const name of options.columns ?? properties.map((p) => p.name)) {
    const prop = properties.find((p) => p.name === name);
    // a bulk value is not a row: a chosen field left empty CLEARS the column, and whether the
    // column tolerates that is the server's per-row answer — nothing here is "required"
    if (prop !== undefined && writable.field(name) === 'editable')
      props.push({...prop, nullable: true});
  }
  if (props.length === 0) {
    notify.warning(`There is no column of ${source.table} you may write.`);
    return null;
  }
  const targets = _targets(source);
  if (targets.length === 0) {
    notify.warning(source.search.peek() === '' ?
      `Select the ${_plural(source)} to write, or filter the list — a bulk edit needs a target.` :
      'A bulk edit matches on the filter alone — clear the search box, or select the rows to write.');
    return null;
  }
  if (!await confirmDiscard(source.session, {action: 'edit rows in bulk'}))
    return null;
  return new Promise((resolve) => _open(source, props, targets, resolve));
}

/** The table's plural name, as the user reads it — never "rows". */
function _plural(source: DomainSource): string {
  return DomainTable.noun(source.schema.info);
}

function _count(source: DomainSource, n: number): string {
  return plural(n, DomainTable.noun(source.schema.info, true), _plural(source));
}

/** The targets a bulk edit can actually write: the frame's selection, and — only when the list is
 * filtered and nothing is searched, which `updateWhere` cannot take — everything it matches. */
function _targets(source: DomainSource): BulkTarget[] {
  const selected = source.selection.peek().length;
  const total = source.total.peek() ?? 0;
  const matching = source.search.peek() === '' && _filterOf(source) !== undefined;
  const targets: BulkTarget[] = [];
  if (selected > 0) {
    const words = DomainTable.noun(source.schema.info, selected === 1);
    targets.push({kind: SELECTED, label: `${selected.toLocaleString()} selected ${words}`, count: selected});
  }
  if (matching)
    targets.push({kind: MATCHING, label: `All ${_count(source, total)} matching the filter`, count: total});
  return targets;
}

function _open(source: DomainSource, props: IProperty[], targets: BulkTarget[],
  resolve: (n: number | null) => void): void {
  const handle = DomainTable.of(source)!;
  const dialog = Dialog.create(`Bulk edit ${_plural(source)}`, {name: 'bulk-edit'});
  dialog.root.classList.add('u2-domain-bulk');
  const form = dialog.runInScope(() =>
    propertyForm(props, {}, {access: source.access.peek().columnPolicy(), layout: 'normal'}));

  // "Apply to" and the field picker head the same form, in its label column
  const header = document.createElement('div');
  header.className = 'u2-domain-bulk-header';
  form.root.prepend(header);
  const target = _applyTo(dialog, form, header, targets);
  const picker = dialog.runInScope(() => new MultiSelect({label: 'Change', name: 'fields', items: [],
    placeholder: 'Choose the fields to change', fill: true}));
  form.add(picker, header);
  // one value on many rows can only collide on a business key
  const keys = new Set(source.schema.info.businessKey);
  dialog.effect(() => picker.setItems(props
    .filter((p) => target.value.count <= 1 || !keys.has(p.name!))
    .map((p) => ({value: p.name!, label: DomainTable.captionOf(p)}))));

  const selection = source.selection.peek();
  for (const prop of props) {
    const input = form.input(prop.name!);
    if (input === undefined)
      continue;
    const current = _current(selection, prop.name!);
    const hint = span('', 'u2-domain-bulk-current');
    input.root.append(hint);
    input.effect(() => {
      input.visible = picker.value.value.includes(prop.name!);
      hint.textContent = Input.isEmpty(input.value.value) ? 'Will be cleared' :
        target.value.kind === SELECTED ? current : '';
    });
  }

  const checks = handle.validators.columns.length > 0;
  const problem = computed(() => {
    const {kind, count} = target.value;
    if (picker.value.value.length === 0)
      return 'Choose the fields to change';
    if (kind === SELECTED && count > MAX_ROWS) {
      return `A bulk edit writes at most ${MAX_ROWS} ${_plural(source)} at a time; ${count} are ` +
        'selected. Narrow the selection, or filter the list and apply to everything it matches.';
    }
    if (kind === MATCHING && count > MAX_ROWS && checks)
      return _tooMany(source);
    return null;
  });
  const caption = computed(() => {
    const {count} = target.value;
    return count > MAX_ROWS ? `Update ${MAX_ROWS} of ${_count(source, count)}` : `Update ${_count(source, count)}`;
  });

  let done = false;
  const finish = (n: number | null) => {
    if (done)
      return;
    done = true;
    resolve(n);
  };
  // busy until the write ends; a refused one keeps what the user chose and typed, the reason in the
  // dialog's error line
  const run = async () => {
    const columns = props.map((p) => p.name!).filter((name) => picker.value.peek().includes(name));
    finish(await _run(source, form, columns, target.peek().kind === SELECTED));
  };
  dialog.effect(() => {
    dialog.onOK(run, {text: caption.value});
  });
  dialog.add(form).okDisabled(problem).onCancel(() => finish(null)).show({modal: true, width: 480});
}

/** "Apply to": plain text where there is one target, a choice where there are two. */
function _applyTo(dialog: Dialog, form: ObjectForm, header: HTMLElement, targets: BulkTarget[]):
  ReadonlySignal<BulkTarget> {
  if (targets.length === 1) {
    header.append(ObjectForm.readonlyField('Apply to', 'target', targets[0].label).row);
    return signal(targets[0]);
  }
  const radio = dialog.runInScope(() => new RadioInput({label: 'Apply to', name: 'target',
    items: targets.map((t) => t.label), value: targets[0].label, nullable: false}));
  form.add(radio, header);
  return computed(() => targets.find((t) => t.label === radio.value.value) ?? targets[0]);
}

/** What a field holds across the selected rows now: their one shared value, or how many differ. */
function _current(rows: readonly RowView[], column: string): string {
  const seen = new Set(rows.map((row) => text(row[Rows.caption(column)] ?? row[column])));
  if (seen.size > 1)
    return `${seen.size} different values`;
  const [only] = seen;
  return only === undefined || only === '' ? 'Currently empty' : `Currently ${only}`;
}

function _tooMany(source: DomainSource): string {
  return `The filter matches more than ${MAX_ROWS} ${_plural(source)}, and a bulk edit checks the ` +
    `table's rules on at most ${MAX_ROWS}. Narrow the filter.`;
}

/** The list's own filter as the seam takes it — `DomainSource` sends exactly this with every
 * query; undefined where nothing is filtered, which the server refuses for an update. */
function _filterOf(source: DomainSource): DomainQueryLike['filter'] {
  const q = source.query.peek();
  return typeof q === 'string' ? (q === '' ? undefined : q) : Filters.toDomainTree(q, {schema: source.schema});
}

async function _run(source: DomainSource, form: ObjectForm, columns: string[], selection: boolean):
  Promise<number> {
  const handle = DomainTable.of(source)!;
  const filter = selection ? Filters.toDomainTree(Filters.group('and',
    [Filters.cond('id', 'in', source.selection.peek().map((row) => row.id))]), {schema: source.schema}) :
    _filterOf(source);
  // a chosen column left empty clears the column; one not chosen is not in the request
  const written = Object.fromEntries(columns.map((name) => {
    const value = form.input(name)?.value.peek();
    return [name, Input.isEmpty(value) ? null : value];
  }));
  const refused = await _check(source, handle, filter, written, selection);
  if (refused !== null)
    throw new Error(refused);
  const report = await handle.table.updateWhere!(filter, written, {limit: MAX_ROWS}).catch((e) => {
    throw new Error(DomainErrors.message(e));
  });
  const what = DomainTable.noun(source.schema.info, report.updated === 1);
  notify.info(`Updated ${report.updated} ${what}${report.hasMore ?
    ` — more match than one edit writes (${MAX_ROWS}); run it again.` : '.'}`);
  await source.refresh();
  return report.updated;
}

/** The table's validators over the target rows as they would land; null where they pass or the
 * table has none. The matching rows are read for it, and refused past what one edit writes. */
async function _check(source: DomainSource, handle: DomainTable, filter: DomainQueryLike['filter'],
  written: Record<string, unknown>, selection: boolean): Promise<string | null> {
  if (handle.validators.columns.length === 0)
    return null;
  const rows = selection ? source.selection.peek() :
    await handle.table.query({filter, limit: MAX_ROWS + 1}) as RowView[];
  if (rows.length > MAX_ROWS)
    return _tooMany(source);
  return firstProblem(handle, rows.map((row) => ({...row, ...written})))?.message ?? null;
}

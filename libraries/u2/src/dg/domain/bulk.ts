/* `domains.bulkEdit` — one value written into many rows at once (WO 3-8), through the platform's
   `POST …/{table}/update`. The dialog offers every column the caller may write, each with an
   include checkbox: only the checked ones reach the server, so "clear this field" (a checked
   column left empty) and "leave it alone" (unchecked) are different requests. The target is the
   selection, or everything the list's filter matches; the server narrows both by the Edit
   predicate and answers how many rows it touched.

   Only targets that can actually be written are offered, and whatever blocks OK is a line above
   the buttons, never a balloon across the screen — a balloon says what landed, nothing else.

   A bulk edit is NOT part of the unit of work — it lands the moment OK is pressed — so the
   session's pending changes are settled through `confirmDiscard` first. */
import {signal, computed} from '../../core/signals.js';
import type {Signal} from '../../core/signals.js';
import {Filters} from '../../core/filter/index.js';
import {div, span} from '../../core/elements.js';
import type {IProperty} from '../../core/property-like.js';
import {Dialog} from '../../components/containers/dialog.js';
import {ChoiceInput} from '../../components/inputs/choice-input.js';
import {notify} from '../../components/display/notify.js';
import {confirmDiscard} from '../../sources/session.js';
import type {DomainSource} from '../../sources/domain-source.js';
import type {DomainQueryLike} from '../../sources/domain-backend.js';
import {propertyForm} from '../forms/object-form.js';
import {DomainTable} from './index.js';
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

/** Writes the checked columns into the rows the dialog's target names; resolves to how many rows
 * the server touched, or null where the dialog was cancelled or the edit refused. */
export async function bulkEdit(source: DomainSource, options: DomainBulkEditOptions = {}): Promise<number | null> {
  const table = DomainTable.of(source)?.table;
  if (table?.updateWhere === undefined) {
    notify.error(`${source.table}: the backend does not support bulk edits`);
    return null;
  }
  // before the narrowing below, which lifts an upper bound rather than reading past a false
  // negative: a trash source is read-only until its rows are restored
  if (source.readOnly.peek()) {
    notify.warning(`${source.table}: deleted rows are read-only until they are restored`);
    return null;
  }
  const access = source.access.peek();
  // column security alone decides what is offered: on a row-mode table the table-level `edit` is a
  // false negative (GOAL "Access"), and the server narrows the update by the row predicate anyway
  const writable = access.narrow({edit: true});
  const properties = source.schema.properties;
  const props: IProperty[] = [];
  for (const name of options.columns ?? properties.map((p) => p.name)) {
    const prop = properties.find((p) => p.name === name);
    // a bulk value is not a row: a checked field left empty CLEARS the column, and whether the
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
  return (source.schema.info.pluralName || 'rows').replace(/_/g, ' ').toLowerCase();
}

/** The targets a bulk edit can actually write: the frame's selection, and — only when the list is
 * filtered and nothing is searched, which `updateWhere` cannot take — everything it matches. */
function _targets(source: DomainSource): {value: string, label: string}[] {
  const selected = source.selection.peek().length;
  const matching = source.search.peek() === '' && _filterOf(source) !== undefined;
  return [
    ...(selected === 0 ? [] : [{value: SELECTED, label: `${selected} selected`}]),
    ...(matching ? [{value: MATCHING, label: `all ${source.total.peek() ?? 0} matching`}] : []),
  ];
}

function _open(source: DomainSource, props: IProperty[], items: {value: string, label: string}[],
  resolve: (n: number | null) => void): void {
  const dialog = Dialog.create(`Bulk edit ${_plural(source)}`, {name: 'bulk-edit'});
  dialog.root.classList.add('u2-domain-bulk');
  const target = dialog.runInScope(() => new ChoiceInput({label: 'Apply to', items,
    value: items[0].value, name: 'target', nullable: false}));
  const values: Record<string, unknown> = {};
  const include = new Map<string, Signal<boolean>>();
  const form = dialog.runInScope(() =>
    propertyForm(props, values, {access: source.access.peek().narrow({edit: true}), layout: 'tall'}));
  for (const prop of props) {
    const on = signal(false);
    include.set(prop.name!, on);
    const input = form.input(prop.name!);
    if (input === undefined)
      continue;
    // OUTSIDE the input's root: a disabled input kills pointer events on everything inside it and
    // disables the controls it holds, so the box that would enable it would be dead in the water
    const row = div([_checkbox(prop, on)], 'u2-domain-bulk-row');
    input.root.replaceWith(row);
    row.append(input.root);
    input.enabled = false;
    input.effect(() => input.enabled = on.value);
  }

  const problem = computed(() => {
    if (!props.some((p) => include.get(p.name!)!.value))
      return 'Check the fields to write.';
    const selected = source.selection.value.length;
    if (target.value.value === SELECTED && selected > MAX_ROWS) {
      return `A bulk edit writes at most ${MAX_ROWS} ${_plural(source)} at a time; ${selected} are ` +
        'selected. Narrow the selection, or filter the list and apply to everything it matches.';
    }
    return null;
  });
  const reason = span('', 'u2-domain-bulk-reason');
  dialog.effect(() => reason.textContent = problem.value ?? '');
  dialog.add(target)
    .add(span('Checked fields are written to every target row; a checked field left empty clears it.',
      'u2-domain-bulk-hint'))
    .add(form).add(reason);

  let done = false;
  const finish = (n: number | null) => {
    if (done)
      return;
    done = true;
    dialog.dispose();
    resolve(n);
  };
  dialog.onOK(() => {
    const columns = props.map((p) => p.name!).filter((name) => include.get(name)!.peek());
    void _run(source, columns, values, target.value.peek() === SELECTED).then(finish);
  }).okEnabled(computed(() => problem.value === null))
    .onCancel(() => finish(null)).show({modal: true, width: 460});
}

/** The include checkbox of one column, in front of the input's whole row: an unchecked column is
 * not in the request at all. */
function _checkbox(prop: IProperty, on: Signal<boolean>): HTMLElement {
  const name = prop.friendlyName ?? prop.name ?? '';
  const caption = `${name.charAt(0).toUpperCase()}${name.slice(1)}`;
  const box = document.createElement('input');
  box.type = 'checkbox';
  box.className = 'u2-input-checkbox u2-domain-bulk-include';
  box.setAttribute('aria-label', caption);
  box.title = `${caption} (leave empty to clear)`;
  box.dataset.u2Include = prop.name;
  box.addEventListener('change', () => on.value = box.checked);
  return box;
}

/** The list's own filter as the seam takes it — `DomainSource` sends exactly this with every
 * query; undefined where nothing is filtered, which the server refuses for an update. */
function _filterOf(source: DomainSource): DomainQueryLike['filter'] {
  const q = source.query.peek();
  return typeof q === 'string' ? (q === '' ? undefined : q) : Filters.toDomainTree(q, {schema: source.schema});
}

async function _run(source: DomainSource, columns: string[], values: Record<string, unknown>,
  selection: boolean): Promise<number | null> {
  const table = DomainTable.of(source)!.table;
  const filter = selection ? Filters.toDomainTree(Filters.group('and',
    [Filters.cond('id', 'in', source.selection.peek().map((row) => row.id))]), {schema: source.schema}) :
    _filterOf(source);
  // a checked column with no value clears the column; an unchecked one is not in the request
  const written = Object.fromEntries(columns.map((name) => [name, values[name] ?? null]));
  try {
    const report = await table.updateWhere!(filter, written, {limit: MAX_ROWS});
    const what = report.updated === 1 ? source.schema.info.singularName.toLowerCase() || 'row' : _plural(source);
    notify.info(`Updated ${report.updated} ${what}${report.hasMore ?
      ` — more match than one edit writes (${MAX_ROWS}); run it again.` : '.'}`);
    await source.refresh();
    return report.updated;
  } catch (e) {
    notify.error(DomainErrors.message(e));
    return null;
  }
}

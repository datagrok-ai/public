/** Aggregation-list editor for `Flow:aggregate`.
 *
 *  "Which columns, aggregated how" is a list, and a list is the one shape a
 *  function signature cannot express with primitives — which is exactly why
 *  `core:Aggregate` takes `List<GroupAggregation>` and ends up unusable on a
 *  canvas. Here it is a JSON string edited as rows of two combos: a column and
 *  an aggregation, both picked, never typed.
 *
 *  The column choices are the upstream table's, filtered by what the chosen
 *  aggregation can actually consume (an average of a string column is a runtime
 *  error waiting to happen), and `count` drops the column combo entirely
 *  because it counts rows.
 *
 *  With no table available the editor does NOT render a row of empty combos —
 *  it says what is missing and, when the upstream is merely uncomputed, offers
 *  to run the slice; the stored aggregations stay visible as text so a flow
 *  loaded from disk is still readable. Same contract as the MPO column mapper,
 *  for the same reason. */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {CustomEditorContext, CustomInputEditor} from '../../utils/func-input-overrides';
import {getParamDisplayName} from '../../utils/dart-proxy-utils';
import {AGGREGATION_TYPES, AggregationSpec, parseAggregations} from '../../ops/data-ops';

/** Aggregations defined only over numbers — offering a string column for them
 *  produces nulls at best. */
const NUMERIC_AGGREGATIONS = new Set([
  'min', 'max', 'sum', 'avg', 'med', 'stdev', 'variance', 'q1', 'q2', 'q3', 'skew', 'kurt',
]);

/** Aggregations that take no column at all. */
const COLUMNLESS_AGGREGATIONS = new Set(['count']);

/** The columns a given aggregation can consume. Pure — unit-tested. */
export function columnsForAggregation(
  type: string, columns: readonly DG.Column[] | null,
): string[] {
  if (!columns || COLUMNLESS_AGGREGATIONS.has(type)) return [];
  const numericOnly = NUMERIC_AGGREGATIONS.has(type);
  return columns.filter((c) => {
    if (!numericOnly) return true;
    try {
      return c.matches('numerical');
    } catch {
      return false;
    }
  }).map((c) => c.name);
}

/** Drop entries whose column is no longer in the table, so switching the
 *  upstream table doesn't leave an aggregation over a column that is gone.
 *
 *  A BLANK column survives: that is a row the user has just added and not
 *  finished, not a stale reference — pruning it would make the Add button
 *  appear to do nothing. The readiness check is what refuses to run on one.
 *  Column-less aggregations always survive. Pure. */
export function pruneAggregations(
  specs: readonly AggregationSpec[], columnNames: readonly string[],
): AggregationSpec[] {
  const known = new Set(columnNames);
  return specs.filter((a) =>
    COLUMNLESS_AGGREGATIONS.has(a.type) || a.column === '' || known.has(a.column));
}

/** The stored form. Blank when there is nothing to store, so an untouched
 *  parameter reads as empty rather than as `[]`. Pure. */
export function serializeAggregations(specs: readonly AggregationSpec[]): string {
  const cleaned = specs
    .filter((a) => a.type !== '')
    .map((a) => ({
      column: COLUMNLESS_AGGREGATIONS.has(a.type) ? '' : a.column,
      type: a.type,
      ...(a.name && a.name.trim() !== '' ? {name: a.name.trim()} : {}),
    }));
  return cleaned.length === 0 ? '' : JSON.stringify(cleaned);
}

/** One-line human form of an aggregation, for the read-only fallback. */
export function describeAggregation(spec: AggregationSpec): string {
  const base = COLUMNLESS_AGGREGATIONS.has(spec.type) ? spec.type : `${spec.type}(${spec.column})`;
  return spec.name && spec.name.trim() !== '' ? `${base} → ${spec.name.trim()}` : base;
}

export function aggregationEditor(param: DG.Property, ctx: CustomEditorContext): CustomInputEditor {
  const ed: CustomInputEditor = {} as CustomInputEditor;
  const host = ui.div([], 'ff-aggregation-editor');
  // The editor replaces the parameter's whole row, so it has to carry the
  // parameter's own name — without it the combos sit under the previous row's
  // label and read as part of it.
  const label = getParamDisplayName(param);
  let specs: AggregationSpec[] = [];

  const commit = (): void => ed.onChanged?.(serializeAggregations(specs));

  /** No columns to choose from — say why, keep what is stored readable. */
  const renderBlocked = (): void => {
    ui.empty(host);
    host.setAttribute('data-blocked', 'true');
    const connected = ctx.isConnected?.('table') ?? false;
    const rows: HTMLElement[] = [ui.divText(label, 'ff-aggregation-editor-label'), ui.divText(connected ?
      'The upstream table has not been computed yet — its columns are what you aggregate.' :
      'Connect a table to choose columns and aggregations.', 'ff-aggregation-editor-blocked-text')];
    if (connected && ctx.produceTable) {
      const load = ui.link('Load it', () => {
        void ctx.produceTable!('table')
          .then((produced) => produced ? render() : grok.shell.error('The flow ran but produced no table.'))
          .catch((e) => grok.shell.error(`Could not load the table: ${e?.message ?? e}`));
      }, 'Run the flow up to the upstream node');
      rows.push(load);
    }
    if (specs.length > 0) {
      rows.push(ui.divText(`Stored: ${specs.map(describeAggregation).join(', ')}`,
        'ff-aggregation-editor-stored'));
    }
    host.appendChild(ui.divV(rows, 'ff-aggregation-editor-blocked'));
  };

  const render = (): void => {
    const columns = ctx.columns('table');
    if (!columns || columns.length === 0) return renderBlocked();

    const names = columns.map((c) => c.name);
    const pruned = pruneAggregations(specs, names);
    if (serializeAggregations(pruned) !== serializeAggregations(specs)) {
      specs = pruned;
      commit();
    } else
      specs = pruned;

    ui.empty(host);
    host.removeAttribute('data-blocked');
    const list = ui.div([], 'ff-aggregation-rows');

    specs.forEach((spec, index) => {
      const columnChoices = columnsForAggregation(spec.type, columns);
      const typeInput = ui.input.choice('', {
        value: spec.type,
        items: [...AGGREGATION_TYPES],
        nullable: false,
        tooltipText: 'How to aggregate',
        onValueChanged: (v) => {
          spec.type = String(v ?? '');
          // The new aggregation may not accept the column that was picked.
          const allowed = columnsForAggregation(spec.type, columns);
          if (spec.column !== '' && !allowed.includes(spec.column)) spec.column = '';
          commit();
          render();
        },
      });
      const rowItems: HTMLElement[] = [typeInput.root];

      if (!COLUMNLESS_AGGREGATIONS.has(spec.type)) {
        const columnInput = ui.input.choice('', {
          value: columnChoices.includes(spec.column) ? spec.column : '',
          items: ['', ...columnChoices],
          nullable: true,
          tooltipText: NUMERIC_AGGREGATIONS.has(spec.type) ?
            'Numeric column to aggregate' : 'Column to aggregate',
          onValueChanged: (v) => {
            spec.column = String(v ?? '');
            commit();
          },
        });
        rowItems.push(columnInput.root);
      }

      const remove = ui.iconFA('times', () => {
        specs.splice(index, 1);
        commit();
        render();
      }, 'Remove this aggregation');
      remove.classList.add('ff-aggregation-remove');
      rowItems.push(remove);

      const row = ui.divH(rowItems, 'ff-aggregation-row');
      row.setAttribute('data-aggregation-index', String(index));
      list.appendChild(row);
    });

    const add = ui.link('Add aggregation', () => {
      // Seed with the first aggregation that has something to aggregate, so a
      // fresh row is immediately meaningful rather than half-blank.
      specs.push({column: '', type: 'avg'});
      commit();
      render();
    }, 'Aggregate one more column');
    add.classList.add('ff-aggregation-add');

    host.appendChild(ui.divV([ui.divText(label, 'ff-aggregation-editor-label'), list, add]));
  };

  // The group-by columns don't change what can be aggregated, but the table
  // input does — and the panel doesn't re-render while focus is inside it.
  ctx.watch('groupByColumns', () => render());

  ed.element = host;
  ed.getValue = (): unknown => serializeAggregations(specs);
  ed.setValue = (v): void => {
    specs = parseAggregations(v);
    render();
  };
  return ed;
}

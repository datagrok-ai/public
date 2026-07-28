/**
 * Compile-time verification of every fenced example in SKILL.md.
 *
 * Every `datagrok-exec` block referenced from SKILL.md appears here as a named
 * function with the same code. If any example uses a non-existent method or
 * wrong arg shape, this file fails to type-check.
 *
 * Each example is wrapped in a function with the same global signature the
 * runtime injects into a `datagrok-exec` block: `(grok, ui, DG, view, t)`.
 * The runtime does NOT inject `grokky` — these examples use the raw js-api.
 */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

type ExampleFn = (
  grok: typeof import('datagrok-api/grok'),
  ui: typeof import('datagrok-api/ui'),
  DG: typeof import('datagrok-api/dg'),
  view: DG.TableView,
  t: DG.DataFrame,
) => Promise<unknown> | unknown;

export const examples: Array<{name: string; fn: ExampleFn}> = [
  {name: 'read-counts', fn: async (_grok, _ui, _DG, _view, t) => {
    return {selected: t.selection.trueCount, unselected: t.selection.falseCount, total: t.rowCount};
  }},

  {name: 'read-selected-indexes', fn: async (_grok, _ui, _DG, _view, t) => {
    return Array.from(t.selection.getSelectedIndexes());
  }},

  {name: 'read-is-row-selected', fn: async (_grok, _ui, _DG, _view, t) => {
    return t.selection.get(17);
  }},

  {name: 'select-predicate-age', fn: async (_grok, _ui, _DG, _view, t) => {
    const age = t.getCol('age');
    t.selection.init((i) => (age.get(i) as number) > 40);
    return {selected: t.selection.trueCount};
  }},

  {name: 'add-predicate-activity', fn: async (_grok, _ui, DG, _view, t) => {
    const act = t.getCol('activity');
    t.selection.or(DG.BitSet.create(t.rowCount, (i) => (act.get(i) as number) > 7));
  }},

  {name: 'intersect-predicate-category', fn: async (_grok, _ui, DG, _view, t) => {
    const cat = t.getCol('category');
    t.selection.and(DG.BitSet.create(t.rowCount, (i) => cat.get(i) === 'A'));
  }},

  {name: 'select-indices-replace', fn: async (_grok, _ui, _DG, _view, t) => {
    const idx = [0, 3, 5, 7, 11];
    t.selection.setAll(false, false);
    for (const i of idx) t.selection.set(i, true, false);
    t.selection.fireChanged();
  }},

  {name: 'remove-indices', fn: async (_grok, _ui, _DG, _view, t) => {
    for (const i of [5, 7]) t.selection.set(i, false, false);
    t.selection.fireChanged();
  }},

  {name: 'clear-selection', fn: async (_grok, _ui, _DG, _view, t) => {
    t.selection.setAll(false);
  }},

  {name: 'select-all', fn: async (_grok, _ui, _DG, _view, t) => {
    t.selection.setAll(true);
  }},

  {name: 'invert-selection', fn: async (_grok, _ui, _DG, _view, t) => {
    t.selection.invert();
  }},

  {name: 'selected-df', fn: async (_grok, _ui, _DG, _view, t) => {
    const subset = t.clone(t.selection);
    return subset;
  }},

  {name: 'selected-df-as-new-view', fn: async (grok, _ui, _DG, _view, t) => {
    const subset = t.clone(t.selection);
    grok.shell.addTableView(subset);
  }},

  {name: 'filter-from-selection', fn: async (_grok, _ui, _DG, _view, t) => {
    t.filter.copyFrom(t.selection);
    return {visible: t.filter.trueCount, selected: t.selection.trueCount};
  }},

  {name: 'selection-from-filter', fn: async (_grok, _ui, _DG, _view, t) => {
    t.selection.copyFrom(t.filter);
    return {selected: t.selection.trueCount};
  }},

  {name: 'set-current-row', fn: async (_grok, _ui, _DG, _view, t) => {
    if (5 >= 0 && 5 < t.rowCount)
      t.currentRowIdx = 5;
  }},

  {name: 'describe-selection', fn: async (_grok, _ui, _DG, _view, t) => {
    const sel = t.selection;
    const all = sel.getSelectedIndexes();
    const MAX_INDEXES = 50;
    const MAX_SAMPLE_COLS = 12;
    const indexes = Array.from(all.slice(0, MAX_INDEXES));

    const out: {
      count: number; total: number; indexes: number[]; currentRowIdx: number;
      sample?: Record<string, unknown>;
    } = {
      count: sel.trueCount,
      total: t.rowCount,
      indexes,
      currentRowIdx: t.currentRowIdx,
    };

    if (all.length > 0) {
      const firstIdx = all[0];
      const cols = t.columns;
      const sample: Record<string, unknown> = {};
      const colCap = Math.min(cols.length, MAX_SAMPLE_COLS);
      for (let c = 0; c < colCap; c++) {
        const col = cols.byIndex(c);
        sample[col.name] = col.get(firstIdx);
      }
      out.sample = sample;
    }
    return out;
  }},

  {name: 'pattern-predicate-most-common', fn: async (_grok, _ui, _DG, _view, t) => {
    const act = t.getCol('activity');
    t.selection.init((i) => (act.get(i) as number) > 7);
  }},

  {name: 'pattern-toggle-single-row', fn: async (_grok, _ui, _DG, _view, t) => {
    const i = 5;
    t.selection.set(i, !t.selection.get(i));
  }},

  {name: 'pattern-top-n-by-column', fn: async (_grok, _ui, _DG, _view, t) => {
    const idx = Array.from(t.getSortedOrder(['activity'], [false])).slice(0, 5);
    t.selection.setAll(false, false);
    for (const i of idx) t.selection.set(i, true, false);
    t.selection.fireChanged();
  }},

  {name: 'pattern-pipeline-filter-select-narrow', fn: async (_grok, _ui, DG, _view, t) => {
    t.selection.copyFrom(t.filter);
    const flag = t.getCol('flag');
    t.selection.and(DG.BitSet.create(t.rowCount, (i) => flag.get(i) === 'ok'));
    t.filter.copyFrom(t.selection);
  }},

  {name: 'subscribe-onSelectionChanged', fn: async (_grok, _ui, _DG, view, t) => {
    const sub = t.onSelectionChanged.subscribe(() => {
      console.log(`selected: ${t.selection.trueCount}`);
    });
    view.subs.push(sub);
  }},

  // Silence "imported but never used" while keeping the imports for type checking.
  {name: '_noop-imports', fn: () => {
    void grok; void ui; void DG;
  }},
];

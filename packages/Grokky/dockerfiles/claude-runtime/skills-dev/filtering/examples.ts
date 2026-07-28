/**
 * Compile-time verification of every fenced example in SKILL.md.
 *
 * Every `datagrok-exec` block referenced from SKILL.md appears here as a
 * named function with the same code. If any example uses a non-existent
 * method or wrong arg shape, this file fails to type-check.
 *
 * Each example is wrapped in a function with the same global signature the
 * runtime injects into a `datagrok-exec` block: `(grok, ui, DG, view, t)`.
 * There is no `grokky` parameter — that wrapper layer was removed; skills
 * teach raw js-api directly.
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
    return {visible: t.filter.trueCount, hidden: t.filter.falseCount, total: t.rowCount};
  }},

  {name: 'read-visible-indexes', fn: async (_grok, _ui, _DG, _view, t) => {
    return t.filter.getSelectedIndexes();
  }},

  {name: 'predicate-keep-age', fn: async (_grok, _ui, _DG, _view, t) => {
    t.filter.init((i) => (t.getCol('age').get(i) as number) > 40);
    return {visible: t.filter.trueCount};
  }},

  {name: 'predicate-multi-column', fn: async (_grok, _ui, _DG, _view, t) => {
    const mw = t.getCol('MW');
    const logP = t.getCol('cLogP');
    t.filter.init((i) => (mw.get(i) as number) < 500 && (logP.get(i) as number) < 5);
  }},

  {name: 'predicate-intersect-existing', fn: async (_grok, _ui, DG, _view, t) => {
    const col = t.getCol('activity');
    const mask = DG.BitSet.create(t.rowCount, (i) => (col.get(i) as number) > 7);
    t.filter.and(mask);
  }},

  {name: 'collaborative-predicate-filter', fn: async (_grok, _ui, DG, view, t) => {
    const sub = t.onRowsFiltering.subscribe((_) => {
      const mw = t.getCol('MW');
      const mask = DG.BitSet.create(t.rowCount, (i) => (mw.get(i) as number) < 500);
      t.filter.and(mask);
    });
    view.subs.push(sub);
    t.rows.requestFilter();
  }},

  {name: 'filter-range-ui', fn: async (_grok, _ui, DG, view, _t) => {
    view.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
      type: DG.FILTER_TYPE.HISTOGRAM, column: 'MW', min: 200, max: 500,
    } as any);
  }},

  {name: 'filter-equals-ui', fn: async (_grok, _ui, DG, view, _t) => {
    view.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
      type: DG.FILTER_TYPE.CATEGORICAL, column: 'category', selected: ['A'],
    } as any);
  }},

  {name: 'filter-in-set-ui', fn: async (_grok, _ui, DG, view, _t) => {
    view.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
      type: DG.FILTER_TYPE.CATEGORICAL, column: 'category', selected: ['A', 'B', 'C'],
    } as any);
  }},

  {name: 'filter-contains-ui', fn: async (_grok, _ui, DG, view, _t) => {
    view.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
      type: DG.FILTER_TYPE.FREE_TEXT, column: 'name', value: 'acid',
    } as any);
  }},

  {name: 'filter-regex-predicate', fn: async (_grok, _ui, _DG, _view, t) => {
    const col = t.getCol('name');
    const re = /acid$/i;
    t.filter.init((i) => {
      const v = col.get(i);
      return v != null && re.test(String(v));
    });
  }},

  {name: 'substructure-smiles-ui', fn: async (_grok, _ui, DG, view, _t) => {
    const query = 'c1ccccc1';
    const molBlock = DG.chem.isMolBlock(query)
      ? query
      : DG.chem.convert(query, DG.chem.isSmarts(query) ? DG.chem.Notation.Smarts : DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
    view.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
      type: DG.FILTER_TYPE.SUBSTRUCTURE,
      column: 'smiles',
      columnName: 'smiles',
      molBlock,
      molBlockFailover: query,
    } as any);
  }},

  {name: 'substructure-smarts-ui', fn: async (_grok, _ui, DG, view, _t) => {
    const query = '[#6][!#1]';
    const molBlock = DG.chem.convert(query, DG.chem.Notation.Smarts, DG.chem.Notation.MolBlock);
    view.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
      type: DG.FILTER_TYPE.SUBSTRUCTURE,
      column: 'smiles',
      columnName: 'smiles',
      molBlock,
      molBlockFailover: query,
    } as any);
  }},

  {name: 'substructure-bitset-fastpath', fn: async (grok, _ui, _DG, _view, t) => {
    const bs = await grok.chem.searchSubstructure(t.col('smiles')!, 'c1ccccc1');
    t.filter.copyFrom(bs);
  }},

  {name: 'combine-via-predicate', fn: async (_grok, _ui, _DG, _view, t) => {
    const mw = t.getCol('MW');
    const logP = t.getCol('cLogP');
    t.filter.init((i) => (mw.get(i) as number) < 500 && (logP.get(i) as number) < 5);
  }},

  {name: 'bitset-chain', fn: async (_grok, _ui, DG, _view, t) => {
    const catA = DG.BitSet.create(t.rowCount, (i) => i % 2 === 0);
    const catB = DG.BitSet.create(t.rowCount, (i) => i % 3 === 0);
    t.filter.setAll(true).and(catA).and(catB);
  }},

  {name: 'clear-filter-view', fn: async (_grok, _ui, _DG, view, t) => {
    t.filter.setAll(true);
    view.getFiltersGroup({createDefaultFilters: false}).setActive(false);
  }},

  {name: 'invert-filter', fn: async (_grok, _ui, _DG, _view, t) => {
    t.filter.invert();
  }},

  {name: 'show-only-hide', fn: async (_grok, _ui, _DG, _view, t) => {
    t.filter.init((i) => t.getCol('category').get(i) === 'X');
    return t.filter.trueCount;
  }},

  {name: 'filtered-df-clone', fn: async (_grok, _ui, _DG, _view, t) => {
    const subset = t.clone(t.filter);
    return subset;
  }},

  {name: 'drop-rows-destructive', fn: async (_grok, _ui, _DG, _view, t) => {
    const before = t.rowCount;
    t.rows.removeWhereIdx((i) => t.getCol('category').get(i) === 'X');
    return {removed: before - t.rowCount};
  }},

  {name: 'drop-hidden-rows', fn: async (_grok, _ui, _DG, _view, t) => {
    t.rows.removeWhereIdx((i) => !t.filter.get(i));
  }},

  {name: 'observe-rows-filtered', fn: async (_grok, _ui, _DG, view, t) => {
    const sub = t.onRowsFiltered.subscribe((_) => console.log(`visible: ${t.filter.trueCount}`));
    view.subs.push(sub);
  }},

  // Silence "imported but never used" while keeping the imports for type checking.
  {name: '_noop-imports', fn: () => {
    void grok; void ui; void DG;
  }},
];

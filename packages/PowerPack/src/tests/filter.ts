import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect, awaitCheck} from '@datagrok-libraries/test/src/test';
import type {FilterBuilderFilter} from '../filter/filter-builder-filter';

category('Filters: builder', () => {
  const query = 'age > 30 and sex = "F"';
  const tree = [{property: 'age', operator: '>', value: 30}, 'and', {property: 'sex', operator: '=', value: 'F'}];
  let df: DG.DataFrame;
  let tv: DG.TableView;
  let fg: DG.FilterGroup;
  let expected: number;

  const filtered = (count: number, timeout: number = 10000) =>
    awaitCheck(() => df.filter.trueCount === count, `expected ${count} rows, got ${df.filter.trueCount}`, timeout);

  // the runner closes every view after a test, so each test opens its own
  function open(): void {
    df = grok.data.demo.demog();
    tv = grok.shell.addTableView(df);
    const age = df.col('age')!;
    const sex = df.col('sex')!;
    expected = 0;
    for (let i = 0; i < df.rowCount; i++) {
      if (age.get(i) > 30 && sex.get(i) === 'F')
        expected++;
    }
  }

  // the panel's first refresh is async (it loads the filter's package) and re-adds every look state
  // that plain add() left there, so add() right after getFiltersGroup() mounts any filter twice —
  // histograms included — hence updateOrAdd; the test bundle carries its own copy of the class, so
  // the filter is found by caption rather than instanceof
  async function addFilter(state: object = {query}): Promise<FilterBuilderFilter> {
    const type = DG.Func.find({package: 'PowerPack', name: 'filterBuilder'})[0].nqName;
    open();
    fg = tv.getFiltersGroup({createDefaultFilters: false});
    fg.updateOrAdd({type, ...state});
    await filtered(expected);
    return fg.filters.find((f) => (f as DG.Filter).caption === 'Filter builder') as FilterBuilderFilter;
  }

  test('query state filters the frame', async () => {
    await addFilter();
    expect(expected > 0 && expected < df.rowCount, true);
    expect(fg.filters.length, 1);
  });

  // the domain wire shape of pre-`{model, query}` layouts is not read back: the filter mounts empty
  test('a tree-only state adds an empty builder', async () => {
    const type = DG.Func.find({package: 'PowerPack', name: 'filterBuilder'})[0].nqName;
    open();
    fg = tv.getFiltersGroup({createDefaultFilters: false});
    fg.updateOrAdd({type, tree});
    await awaitCheck(() => fg.filters.length === 1, 'the filter was not added', 5000);
    expect((fg.filters[0] as FilterBuilderFilter).saveState().query, '');
    expect(df.filter.trueCount, df.rowCount);
  });

  test('saveState carries model and query', async () => {
    const state = (await addFilter()).saveState();
    expect(state.query, query);
    expect('tree' in state, false);
    expect(state.model.op, 'and');
    expect(state.model.nodes.map((n: any) => `${n.property} ${n.operator} ${n.value}`).join(', '), 'age > 30, sex = F');
  });

  test('a model state through the platform filters the frame', async () => {
    const model = {op: 'and', nodes: [{property: 'age', operator: '>', value: 30}, {property: 'sex', operator: '=', value: 'F'}]};
    const f = await addFilter({model});
    expect(f.saveState().query, query);
  });

  test('an incomplete row leaves the complete ones filtering', async () => {
    const f = await addFilter();
    f.builder!.addCondition();
    await awaitCheck(() => f.filterSummary === query && f.isFiltering, 'the complete rows still filter', 5000);
    await filtered(expected);
    expect(f.saveState().query.startsWith(`${query} and `), true);
  });

  test('the status line is off by default; a state turns it on and saves it', async () => {
    const f = await addFilter();
    const status = f.root.querySelector<HTMLElement>('.power-pack-filter-builder-status')!;
    expect(status.hidden, true);
    expect(f.saveState().showStatus, false);
    f.applyState({query, showStatus: true});
    await awaitCheck(() => !status.hidden && status.textContent === query, 'the status line did not show the query', 5000);
    expect(f.saveState().showStatus, true);
  });

  test('applyState restores from a query', async () => {
    const f = await addFilter();
    f.applyState({query: 'age > 200'});
    await filtered(0);
    f.applyState({query});
    await filtered(expected);
    expect(f.saveState().query, query);
  });

  test('a layout round trip keeps the filter', async () => {
    await addFilter();
    const layout = tv.saveLayout();
    open();
    tv.loadLayout(layout);
    await filtered(expected);
    const restored = tv.getFiltersGroup({createDefaultFilters: false});
    expect(restored.filters.length, 1);
    expect((restored.filters[0] as FilterBuilderFilter).saveState().query, query);
  });

  // Chem's `Contains` has no grammar spelling and no domain form: only `model` carries it across a
  // layout save/load; the semType is set by hand — detection is asynchronous and this test does
  // not wait for it
  test('a substructure filter survives a layout round trip', async () => {
    const smiles = 'c1ccccc1';
    const load = async () => {
      df = await grok.data.files.openTable('System:DemoFiles/chem/smiles.csv');
      df.col('canonical_smiles')!.semType = DG.SEMTYPE.MOLECULE;
      tv = grok.shell.addTableView(df);
    };
    await load();
    const out: DG.Column = await grok.functions.call('Chem:searchSubstructure',
      {molStringsColumn: df.col('canonical_smiles'), molString: smiles, molBlockFailover: ''});
    const hits: number = out.get(0).trueCount;
    expect(hits > 0 && hits < df.rowCount, true);
    const type = DG.Func.find({package: 'PowerPack', name: 'filterBuilder'})[0].nqName;
    fg = tv.getFiltersGroup({createDefaultFilters: false});
    fg.updateOrAdd({type, model: {op: 'and', nodes: [{property: 'canonical_smiles', operator: 'Contains', value: smiles}]}});
    await filtered(hits, 30000);
    const state = (fg.filters.find((f) => (f as DG.Filter).caption === 'Filter builder') as FilterBuilderFilter).saveState();
    expect('tree' in state, false);
    expect(state.model.nodes[0].operator, 'Contains');
    const layout = tv.saveLayout();
    await load();
    tv.loadLayout(layout);
    await filtered(hits, 30000);
    const restored = tv.getFiltersGroup({createDefaultFilters: false});
    expect(restored.filters.length, 1);
    expect((restored.filters[0] as FilterBuilderFilter).saveState().model.nodes[0].value, smiles);
  });

  // the sex row's operator picker: what the user is offered for that column
  test('a semantic type set after attach re-derives the operators', async () => {
    const f = await addFilter();
    const sexOps = () => Array.from(f.root.querySelectorAll('[data-u2="filter-row"]')[1]
      .querySelectorAll<HTMLOptionElement>('[data-u2-part="op"] option')).map((o) => o.value);
    expect(sexOps().includes('Contains'), false);
    const age = df.col('age')!;
    let over30 = 0;
    for (let i = 0; i < df.rowCount; i++) {
      if (age.get(i) > 30)
        over30++;
    }
    df.col('sex')!.semType = DG.SEMTYPE.MOLECULE;
    await awaitCheck(() => sexOps().includes('Contains'), 'the sex row did not gain the molecule operator', 5000);
    // the row keeps `=` as its (now inapplicable) pick and stops filtering; the age row goes on
    expect(sexOps().includes('like'), false);
    await filtered(over30);
    expect(f.saveState().query, query);
    df.col('sex')!.semType = '';
    await awaitCheck(() => sexOps().includes('like'), 'the sex row did not get its string operators back', 5000);
    await filtered(expected);
  });

  test('detach clears the filter', async () => {
    const f = await addFilter();
    fg.remove(f);
    await awaitCheck(() => f.isDetached && !f.isFiltering && f.bitset === null, 'not detached', 5000);
    await filtered(df.rowCount);
  });
});

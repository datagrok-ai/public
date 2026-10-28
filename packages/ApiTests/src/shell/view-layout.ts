import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test, waitUntil} from '@datagrok-libraries/test/src/test';

/** A view state with a scatter plot and a Filters viewer over demog, read once the filters reached it. */
async function scatterAndFilters(): Promise<string> {
  const df = grok.data.demo.demog(20);
  const tv = grok.shell.addTableView(df);
  try {
    tv.scatterPlot({x: 'height', y: 'weight', color: 'race'});
    const filters = tv.getFiltersGroup({createDefaultFilters: false});
    filters.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'sex'});
    filters.updateOrAdd({type: DG.FILTER_TYPE.HISTOGRAM, column: 'age'});
    let state = '';
    await waitUntil(() => (state = tv.saveLayout().viewState).includes('"age"'),
      'the filters never reached the layout', 10000);
    return state;
  } finally {
    tv.close();
    grok.shell.closeTable(df);
  }
}

/** A view state with one viewer of this type and these settings, as a stored layout keeps it. */
function viewerState(type: string, look: {[key: string]: unknown}): string {
  return JSON.stringify({containerType: 'fill', children: [{containerType: 'panel', state: {element: {type, look}}}]});
}

category('Shell: ViewLayout', () => {
  test('boundColumnNames of a view state lists the viewers\' and the filters\' columns', async () => {
    const layout = DG.ViewLayout.fromViewState(await scatterAndFilters());
    expectArray([...layout.boundColumnNames].sort(), ['age', 'height', 'race', 'sex', 'weight']);
    expect(layout.columns.length, 0, 'a layout from a view state carries no columns of its own');
  });

  test('boundColumnNames of a stored layout needs no open view or table', async () => {
    const json = DG.ViewLayout.fromViewState(await scatterAndFilters()).toJson();
    expectArray([...DG.ViewLayout.fromJson(json).boundColumnNames].sort(), ['age', 'height', 'race', 'sex', 'weight']);
  });

  test('boundColumnNames names a column once, whatever its case', async () => {
    const state = viewerState('Scatter plot', {xColumnName: 'Height', yColumnName: 'HEIGHT'});
    expectArray(DG.ViewLayout.fromViewState(state).boundColumnNames, ['Height']);
  });

  test('boundColumnNames reads a viewer type nothing registered from its stored settings', async () => {
    const state = viewerState('ApiTests unregistered viewer',
      {xColumnName: 'age', valueColumnNames: ['height', 'weight'], title: 'not a column'});
    expectArray([...DG.ViewLayout.fromViewState(state).boundColumnNames].sort(), ['age', 'height', 'weight']);
  });

  test('a grid binds no column', async () => {
    const df = grok.data.demo.demog(5);
    const tv = grok.shell.addTableView(df);
    try {
      expectArray(DG.ViewLayout.fromViewState(tv.saveLayout().viewState).boundColumnNames, []);
    } finally {
      tv.close();
      grok.shell.closeTable(df);
    }
  });
}, {owner: 'askalkin@datagrok.ai'});

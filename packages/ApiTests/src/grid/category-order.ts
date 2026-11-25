import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test, TestOptions, waitUntil} from '@datagrok-libraries/test/src/test';

const options: TestOptions = {};

category('Viewers: category order', () => {
  before(async () => {
    if (!DG.Viewer.scatterPlot(grok.data.demo.demog(10)).props.hasProperty('xCategoryOrder'))
      options.skipReason = 'the platform has no category order yet';
  });

  test('category order round-trip', async () => {
    const demog = grok.data.demo.demog(100);
    const byAge = (aggr: string) => JSON.stringify({column: 'age', aggr: aggr, order: 'desc'});
    const aggrOf = (viewer: DG.Viewer, name: string) => JSON.parse(viewer.getOptions(true).look[name]).aggr;
    const sp = demog.plot.scatter({x: 'race', y: 'sex'});
    sp.setOptions({xCategoryOrder: byAge('avg'), yCategoryOrder: byAge('max')});
    expect(aggrOf(sp, 'xCategoryOrder'), 'avg');
    expect(aggrOf(sp, 'yCategoryOrder'), 'max');
    const lc = demog.plot.line({x: 'race'});
    lc.setOptions({xCategoryOrder: byAge('min')});
    expect(aggrOf(lc, 'xCategoryOrder'), 'min');
    const bp = demog.plot.box({category1ColumnName: 'race'});
    bp.setOptions({categoryOrder: byAge('sum')});
    expect(aggrOf(bp, 'categoryOrder'), 'sum');
  }, options);

  test('property keys computed for distinct categories', async () => {
    const calls: number[] = [];
    // grok.functions has no unregister; the func stays registered, harmless under the private ApiTestsNames semType
    const f = grok.functions.register({
      signature: 'column apiTestsNameLength(column/ApiTestsNames values, string metric)',
      tags: 'categoryOrderer',
      run: (values: any, _metric: string) => {
        const col: DG.Column = DG.toJs(values);
        calls.push(col.length);
        return DG.Column.fromList('double', 'keys',
          col.toList().map((s: string | null) => (s ?? '').length));
      },
    });
    const df = DG.DataFrame.fromCsv('names,val\naaa,1\nbb,2\nc,3\naaa,4\nbb,5');
    df.col('names')!.semType = 'ApiTestsNames';
    const view = grok.shell.addTableView(df);
    try {
      const sp = view.scatterPlot({x: 'names', y: 'val'});
      sp.setOptions({xCategoryOrder: JSON.stringify({property: `${f.nqName}(length)`})});
      // the keys arrive asynchronously, and the calculator only ever sees the distinct categories
      await waitUntil(() => calls.length > 0, 'the calculator was not called');
      // a recompute (filter, data change) may repeat the call - every one of them sees the categories
      expect(calls.every((n) => n === 3), true, `distinct categories, not rows: [${calls}]`);
    } finally {
      view.close();
      grok.shell.closeTable(df);
    }
  }, options);
});

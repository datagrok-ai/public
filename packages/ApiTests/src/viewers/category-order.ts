import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test, TestOptions, waitUntil} from '@datagrok-libraries/test/src/test';

const nodeSkip = typeof process !== 'undefined' ? 'browser-only: a table view and grok.functions.register' : undefined;
const anywhere: TestOptions = {};
const browserOnly: TestOptions = {skipReason: nodeSkip};

category('Viewers: category order', () => {
  before(async () => {
    if (!DG.Viewer.scatterPlot(grok.data.demo.demog(10)).props.hasProperty('xCategoryOrder')) {
      for (const options of [anywhere, browserOnly])
        options.skipReason ??= 'the platform has no category order yet';
    }
  });

  test('category order round-trip', async () => {
    const demog = grok.data.demo.demog(100);
    const byProperty = (property: string) => JSON.stringify({property: property, aggr: 'avg', order: 'desc'});
    const propertyOf = (look: {[key: string]: any}, name: string) => JSON.parse(look[name]).property;
    const sp = demog.plot.scatter();
    sp.setOptions({xCategoryOrder: byProperty('Chem:molecularProperty(MW)'),
      yCategoryOrder: byProperty('Chem:molecularProperty(LogP)')});
    const spLook = sp.getOptions(true).look;
    expect(propertyOf(spLook, 'xCategoryOrder'), 'Chem:molecularProperty(MW)');
    expect(propertyOf(spLook, 'yCategoryOrder'), 'Chem:molecularProperty(LogP)');
    const lc = demog.plot.line();
    lc.setOptions({xCategoryOrder: byProperty('Chem:molecularProperty(MW)')});
    expect(propertyOf(lc.getOptions(true).look, 'xCategoryOrder'), 'Chem:molecularProperty(MW)');
    const bp = demog.plot.box();
    bp.setOptions({categoryOrder: byProperty('Chem:molecularProperty(MW)')});
    expect(propertyOf(bp.getOptions(true).look, 'categoryOrder'), 'Chem:molecularProperty(MW)');
  }, anywhere);

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
  }, browserOnly);
});

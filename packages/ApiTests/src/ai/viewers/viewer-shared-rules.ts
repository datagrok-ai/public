import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, settled, test} from '@datagrok-libraries/test/src/test';
import {demog, look, until, withTableView} from '../helpers';

// Rules every viewer shares, checked once over the viewers here: the column rename has one code path for all of
// them (DataFrameViewer._renameColumnInLook, JS viewers included). No rows, deleted rows, a removed column,
// non-finite values and layouts naming a missing column are the loops of the Dart viewers tests; out-of-bounds
// values, Property.clamp's.
category('AI: Viewers: Shared rules', () => {
  const renamed: {[name: string]: string} = {race: 'origin', age: 'years'};
  const afterRename = (value: any): any => Array.isArray(value) ? value.map(afterRename) : renamed[value] ?? value;

  function viewerCases(): [string, {[name: string]: any}][] {
    const cases: [string, {[name: string]: any}][] = [
      [DG.VIEWER.SCATTER_PLOT, {xColumnName: 'age', yColumnName: 'height', colorColumnName: 'race'}],
      [DG.VIEWER.PIE_CHART, {categoryColumnName: 'race'}],
      [DG.VIEWER.PC_PLOT, {columnNames: ['age', 'height', 'weight']}],
      [DG.VIEWER.TILE_VIEWER, {lanesColumnName: 'race'}],
      [DG.VIEWER.TREE_MAP, {splitByColumnNames: ['race', 'sex']}],
      [DG.VIEWER.TRELLIS_PLOT, {xColumnNames: ['sex'], yColumnNames: ['race']}],
      [DG.VIEWER.HISTOGRAM, {valueColumnName: 'age', splitColumnName: 'race'}],
      [DG.VIEWER.BAR_CHART, {valueColumnName: 'age', splitColumnName: 'race'}],
      [DG.VIEWER.BOX_PLOT, {valueColumnName: 'age', categoryColumnNames: ['race']}],
      [DG.VIEWER.LINE_CHART, {xColumnName: 'age', yColumnNames: ['height']}],
    ];
    const charts = DG.Func.find({package: 'Charts', meta: {role: DG.FUNC_TYPES.VIEWER}}).map((f) => f.friendlyName);
    if (charts.includes('Sunburst'))
      cases.push(['Sunburst', {hierarchyColumnNames: ['sex', 'race']}]);
    if (charts.includes('Word cloud'))
      cases.push(['Word cloud', {wordColumnName: 'race'}]);
    return cases;
  }

  test('a renamed column is followed by every setting that named it', async () => {
    const df = demog(100);
    await withTableView(df, async (tv) => {
      const viewers = viewerCases().map(([type, options]) => ({type, options, viewer: tv.addViewer(type, options)}));
      // a JS viewer's host reports nothing pending before its viewer has loaded, and its options land then
      for (const v of viewers) {
        await until(() => Object.entries(v.options)
          .every(([name, value]) => JSON.stringify(look(v.viewer)[name]) === JSON.stringify(value)));
        await settled(v.viewer);
      }
      const shown = viewers.map((v) => v.viewer.getWidgetStatus().values?.['rows shown']);
      grok.shell.clearLastError();
      df.col('race')!.name = 'origin';
      df.col('age')!.name = 'years';
      for (const v of viewers)
        await settled(v.viewer);
      for (const [i, v] of viewers.entries()) {
        for (const [name, value] of Object.entries(v.options))
          expect(JSON.stringify(look(v.viewer)[name]), JSON.stringify(afterRename(value)), `${v.type} ${name}`);
        const after = v.viewer.getWidgetStatus().values?.['rows shown'];
        expect(after === shown[i], true, `${v.type} rows shown ${shown[i]} before the rename, ${after} after it`);
      }
      expect(await grok.shell.lastError ?? '', '');
    });
  });

  test('a renamed column is followed in a viewer\'s Filter formula (core#1213)', async () => {
    const df = demog(100);
    await withTableView(df, async (tv) => {
      const plot = tv.addViewer(DG.VIEWER.SCATTER_PLOT,
        {xColumnName: 'age', yColumnName: 'height', filter: '${age} > 44'});
      await settled(plot);
      const shown = plot.getWidgetStatus().values?.['rows shown'];
      df.col('age')!.name = 'years';
      await settled(plot);
      expect(look(plot)['filter'], '${years} > 44');
      expect(plot.getWidgetStatus().values?.['rows shown'], shown);
    });
  });
}, {owner: 'drizhinashvili@datagrok.ai'});

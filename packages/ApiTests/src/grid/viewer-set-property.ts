import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';


category('Viewers: setProperty', () => {
  let demog: DG.DataFrame;
  let barChart: DG.Viewer;

  before(async () => {
    demog = grok.data.demo.demog(100);
    barChart = demog.plot.bar();
  });

  test('string property', async () => {
    barChart.setOptions({valueAggrType: 'avg'});
    const options = barChart.getOptions(true);
    expect(options.look['valueAggrType'], 'avg');
  });

  test('bool property', async () => {
    barChart.setOptions({showValueAxis: false});
    const options = barChart.getOptions(true);
    expect(options.look['showValueAxis'], false);
  });

  test('int property', async () => {
    barChart.setOptions({maxCategoryWidth: 200});
    const options = barChart.getOptions(true);
    expect(options.look['maxCategoryWidth'], 200);
  });

  test('list property (colorScheme)', async () => {
    const scheme = [4294922560, 4294944000, 4283477800];
    barChart.setOptions({linearColorScheme: scheme});
    const options = barChart.getOptions(true);
    expectArray(options.look['linearColorScheme'], scheme);
  });

  test('string property does not call mapToObject', async () => {
    barChart.setOptions({barSortOrder: 'desc'});
    const options = barChart.getOptions(true);
    expect(options.look['barSortOrder'], 'desc');
  });

  test('bool property does not call mapToObject', async () => {
    barChart.setOptions({relativeValues: true});
    const options = barChart.getOptions(true);
    expect(options.look['relativeValues'], true);
  });
});

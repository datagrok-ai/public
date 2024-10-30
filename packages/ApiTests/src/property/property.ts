import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, test, expect, expectExceptionAsync, expectArray, before} from '@datagrok-libraries/utils/src/test';

category('Property: General', () => {
  test('eval', async () => {
    const dfList: DG.DataFrame[] = await grok.functions
      .eval('OpenServerFile("System:AppData/ApiTests/datasets/demog.csv")');
    expect(dfList[0].columns instanceof DG.ColumnList, true);
  });

  test('call', async () => {
    const dfList: DG.DataFrame[] = await grok.functions
      .call('OpenServerFile', {'fullPath': 'System:AppData/ApiTests/datasets/demog.csv'});
    expect(dfList[0].columns instanceof DG.ColumnList, true);
  });
  
  test('def param', async () => {
    await grok.functions.call('AddNewColumn', {table: grok.data.demo.demog(), expression: 'test', name: 'test'});
  });

  test('call params', async () => {
    expect(await grok.functions.call('sin', {y: 0.5}), null, '{y: 0.5}');
    expect(await grok.functions.call('sin', {x: 0.5}), 0.479425538604203, '{x: 0.5}');
    expect(await grok.functions.call('sin', {X: 0.5}), 0.479425538604203, '{X: 0.5}');
    await expectExceptionAsync(() => grok.functions.call('qqqqqq'));
  });

  test('query params', async () => {
    // expect(await grok.data.query('ApiTests:dummyPackageQuery', {y: 0.5}), null, '{y: 0.5}');
    expect((await grok.data.query('ApiTests:dummyPackageQuery', {x: 0.5})).get('res', 0), 0.5, '{x: 0.5}');
    expect((await grok.data.query('ApiTests:dummyPackageQuery', {X: 0.5})).get('res', 0), 0.5, '{X: 0.5}');
    await expectExceptionAsync(() => grok.data.query('ApiTests:qqqqqq'));
  });
});

category('Property: Header parsing', () => {
  let func!: DG.Func;
  before(async () => {
    func = await grok.functions.eval(`ApiTests:propertyParsing`);
  });

  test('Choices', async () => {
    expectArray(func.inputs[0].choices, ['Standardized', 'Actual']);
  }, {skipReason: 'https://reddata.atlassian.net/browse/GROK-15701'});
  test('Default string val', async () => {
    expect(JSON.parse(func.inputs[1].options['default']), 'Default val');
  });
  test('Default int val', async () => {
    expect(JSON.parse(func.inputs[2].options['default']), 3);
  });
  test('Caption', async () => {
    expect(func.inputs[3].caption, 'My custom caption');
  });
  test('Category', async () => {
    expect(func.inputs[4].category, 'My custom category');
  });
  test('showSlider', async () => {
    expect(func.inputs[5].showSlider, true);
  });
  test('showPlusMinus', async () => {
    expect(func.inputs[6].showPlusMinus, true);
  });
  test('step', async () => {
    expect(func.inputs[7].step, 2);
  });
  test('min', async () => {
    expect(func.inputs[8].min, 1);
  });
  test('max', async () => {
    expect(func.inputs[8].max, 10);
  });
  test('format', async () => {
    expect(func.inputs[9].format, '#.000');
  });
  test('description', async () => {
    expect(func.inputs[10].description, 'My best description');
  });
  test('Complex caption', async () => {
    expect(func.inputs[11].caption, 'MIC O2 P/V exponent');
    expect(ui.input.forProperty(func.inputs[11]).caption, 'MIC O2 P/V exponent');
  }, {skipReason: 'https://reddata.atlassian.net/browse/GROK-15768'});
});

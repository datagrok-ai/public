import * as grok from 'datagrok-api/grok';
import {category, expect, expectObject, expectTable, test} from '@datagrok-libraries/utils/src/test';
import {DataFrame, FileSource} from 'datagrok-api/dg';
import {getScriptTime} from './python_functions';
import dayjs from 'dayjs';

category('Scripts: R scripts', async () => {
  test('int, double, bool, string input/output', async () => {
    const int = 2;
    const double = 0.3;
    const bool = true;
    const str = 'Datagrok';
    const result = await grok.functions.call('CVMTests:RSimple',
      {'integer_input': int, 'double_input': double, 'bool_input': bool, 'string_input': str});
    expectObject(result, {'integer_output': int, 'double_output': double,
      'bool_output': bool, 'string_output': str});
  });

  test('Dataframe input/output', async () => {
    const df = DataFrame.fromCsv(`id,date,name\nid1,${Date.now()},datagrok`);
    const result = await grok.functions.call('CVMTests:RDataframe',
      {'df': df, 'dfNumerical': df, 'dfCategorical': df});
    expectTable(result['resultDf'], df);
    expectTable(result['resultNumerical'], df);
    expectTable(result['resultCategorical'], df);
  });

  test('Graphics output, Column input', async () => {
    const df = DataFrame.fromCsv('x,y\n1,2\n3,4\n5,6');
    const result = await grok.functions.call('CVMTests:RGraphics',
      {'df': df, 'xName': 'x', 'yName': 'y'});
    expect(result && result.length > 0, true);
  });

  test('Dataframe performance test 15 consequently', async () => {
    const results = [];
    for (let i = 0; i < 15; i++) {
      results.push(await getScriptTime('CVMTests:RSingleDf',
        {'df': grok.data.demo.demog(10000)}));
    }
    const sum = results.reduce((p, c) => p + c, 0);
    return {'Average time': sum / results.length,
      'Min time': Math.min(...results), 'Max time': Math.max(...results)};
  }, {timeout: 120000});

  test('Dataframe performance test 25 parallel', async () => {
    const calls = [];
    for (let i = 0; i < 25; i++) {
      calls.push(getScriptTime('CVMTests:RSingleDf',
        {'df': grok.data.demo.demog(10000)}));
    }
    const results = await Promise.all(calls);
    const sum = results.reduce((p, c) => p + c, 0);
    return {'Average time': sum / results.length,
      'Min time': Math.min(...results), 'Max time': Math.max(...results)};
  }, {timeout: 120000});

  test('Calculated column test', async () => {
    const df = DataFrame.fromCsv('x,y\n1,2\n3,4\n5,6');
    const column = await df.columns.addNewCalculated('new', 'CVMTests:RCalcColumn(${x} + ${y})');
    expect(df.columns.contains(column.name), true);
    expect(column.get(0), 6);
  });

  test('Calculated column performance', async () => {
    const df = grok.data.demo.demog(10000);
    const start = Date.now();
    await df.columns.addNewCalculated('new', 'CVMTests:RCalcColumn(${age})');
    return `Execution time: ${Date.now() - start}`;
  });

  test('Datetime input/output', async () => {
    const now = dayjs();
    const result = await grok.functions.call('CVMTests:RDate',
      {'input_datetime': now});
    expect(result, now);
  });

  test('File input', async () => {
    const files = await new FileSource('System:AppData/CVMTests').list('datasets', false, 'cars');
    const result = await grok.functions.call('CVMTests:RLinesCount',
      {'file': files[0], 'header': true, 'separator': ',', 'dec': '.'});
    expect(result, 30);
  });

  test('Column list', async () => {
    const df = DataFrame.fromCsv(`id,date,name\nid1,${Date.now()},datagrok`);
    const cols = df.columns.toList();
    cols.splice(2);
    const result = await grok.functions.call('CVMTests:RColumnList',
      {'df': df, 'cols': cols});
    df.columns.remove('name');
    expectTable(result, df);
  });
});

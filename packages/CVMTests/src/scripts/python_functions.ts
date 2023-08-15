import * as grok from 'datagrok-api/grok';
import {category, expect, expectObject, expectTable, test} from '@datagrok-libraries/utils/src/test';
import {DataFrame, FileSource} from 'datagrok-api/dg';
import dayjs from 'dayjs';

category('Scripts: Python scripts', async () => {
  test('int, double, bool, string input/output', async () => {
    const int = 2;
    const double = 0.3;
    const bool = true;
    const str = 'Datagrok';
    const result = await grok.functions.call('CVMTests:PythonSimple',
      {'integer_input': int, 'double_input': double, 'bool_input': bool, 'string_input': str});
    expectObject(result, {'integer_output': int, 'double_output': double,
      'bool_output': bool, 'string_output': str});
  });

  test('int, double, bool, string default values input/output', async () => {
    const result = await grok.functions.call('CVMTests:PythonSimple');
    expectObject(result, {'integer_output': 1, 'double_output': 3.14,
      'bool_output': true, 'string_output': '4.14string input'});
  }, {skipReason: 'Default parameters are not working in JS-API'});

  test('String choices default', async () => {
    const result = await grok.functions.call('CVMTests:PythonStringChoices');
    expect(result, 'int');
  }, {skipReason: 'Default parameters are not working in JS-API'});

  test('String choices valid', async () => {
    const result = await grok.functions.call('CVMTests:PythonStringChoices',
      {'choices': 'boolean'});
    expect(result, 'boolean');
  });

  test('Dataframe input/output', async () => {
    const df = DataFrame.fromCsv(`id,date,name\nid1,${Date.now()},datagrok`);
    const result = await grok.functions.call('CVMTests:PythonDataframe',
      {'df': df, 'dfNumerical': df, 'dfCategorical': df});
    expectTable(result['resultDf'], df);
    expectTable(result['resultNumerical'], df);
    expectTable(result['resultCategorical'], df);
  });

  test('Graphics output, Column input', async () => {
    const df = DataFrame.fromCsv('x,y\n1,2\n3,4\n5,6');
    const result = await grok.functions.call('CVMTests:PythonGraphics',
      {'df': df, 'xName': 'x', 'yName': 'y'});
    expect(result && result.length > 0, true);
  });

  test('Dataframe performance test 15 consequently', async () => {
    const results = [];
    for (let i = 0; i < 15; i++) {
      results.push(await getScriptTime('CVMTests:PythonSingleDf',
        {'df': grok.data.demo.demog(10000)}));
    }
    const sum = results.reduce((p, c) => p + c, 0);
    return {'Average time': sum / results.length,
      'Min time': Math.min(...results), 'Max time': Math.max(...results)};
  }, {timeout: 120000});

  test('Dataframe performance test 25 parallel', async () => {
    const calls = [];
    for (let i = 0; i < 25; i++) {
      calls.push(getScriptTime('CVMTests:PythonSingleDf',
        {'df': grok.data.demo.demog(10000)}));
    }
    const results = await Promise.all(calls);
    const sum = results.reduce((p, c) => p + c, 0);
    return {'Average time': sum / results.length,
      'Min time': Math.min(...results), 'Max time': Math.max(...results)};
  }, {timeout: 120000});

  test('Calculated column test', async () => {
    const df = DataFrame.fromCsv('x,y\n1,2\n3,4\n5,6');
    const column = await df.columns.addNewCalculated('new', 'CVMTests:PythonCalcColumn(${x} + ${y})');
    expect(df.columns.contains(column.name), true);
    expect(column.get(0), 6);
  });

  test('Calculated column performance', async () => {
    const df = grok.data.demo.demog(10000);
    const start = Date.now();
    await df.columns.addNewCalculated('new', 'CVMTests:PythonCalcColumn(${age})');
    return `Execution time: ${Date.now() - start}`;
  });

  test('Column list', async () => {
    const df = DataFrame.fromCsv(`id,date,name\nid1,${Date.now()},datagrok`);
    const cols = df.columns.toList();
    cols.splice(2);
    const result = await grok.functions.call('CVMTests:PythonColumnList',
      {'df': df, 'cols': cols});
    df.columns.remove('name');
    expectTable(result, df);
  }, {skipReason: 'It doesn\'t work'});

  test('Map type input/output', async () => {
    const result = await grok.functions.call('CVMTests:PythonMap',
      {'input_map': {'hello': 'world'}, 'unique_key': 'my_key'});
    expectObject(result, {'hello': 'world', 'my_key': 'Datagrok'});
  });

  test('Environment string', async () => {
    const result = await grok.functions.call('CVMTests:PythonAnchorsCount',
      {'html': '  <ul>\n' +
              '      <li><a href="https://example.com">Website</a></li>\n' +
              '  <li><a href="mailto:m.bluth@example.com">Email</a></li>\n' +
              '  <li><a href="tel:+123456789">Phone</a></li>\n' +
              '  </ul>'});
    expect(result, 3);
  }, {timeout: 120000});

  test('File type input', async () => {
    const files = await new FileSource('System:AppData/CVMTests').list('datasets', false, 'cars');
    const result = await grok.functions.call('CVMTests:PythonLinesCount',
      {'file': files[0]});
    expect(79498, result);
  }, {skipReason: 'File input is not working'});

  test('File type input and environment yaml', async () => {
    const files = await new FileSource('System:AppData/CVMTests').list('images', false, 'silver.jpg');
    const result = await grok.functions.call('CVMTests:ImagePixelCount',
      {'fileInput': files[0]});
    expect(79498, result);
  }, {timeout: 120000, skipReason: 'File input is not working'});

  test('Datetime input/output', async () => {
    const now = dayjs();
    const result = await grok.functions.call('CVMTests:PythonDate',
      {'input_datetime': now});
    expect(result, now);
  });
});

export async function getScriptTime(name: string, params: object = {}): Promise<number> {
  const start = Date.now();
  await grok.functions.call(name,
    params);
  return Date.now() - start;
}

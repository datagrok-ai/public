import * as grok from 'datagrok-api/grok';
import {category, test, expect, expectObject, expectTable} from '@datagrok-libraries/utils/src/test';
import dayjs from 'dayjs';
import {DataFrame} from 'datagrok-api/dg';

const langs = ['Python', 'R', 'Julia', 'NodeJS', 'Octave', 'Grok', 'JavaScript'];

for (const lang of langs) {
  category(`Scripts: ${lang} scripts`, () => {
    test('int, double, bool, string input/output', async () => {
      const int = 2;
      const double = 0.3;
      const bool = true;
      const str = 'Datagrok';
      const result = await grok.functions.call(`CVMTests:${lang}Simple`,
        {'integer_input': int, 'double_input': double, 'bool_input': bool, 'string_input': str});
      expectObject(result, {'integer_output': int, 'double_output': double,
        'bool_output': bool, 'string_output': str});
    });

    test('Datetime input/output', async () => {
      const now = dayjs();
      const result = await grok.functions.call(`CVMTests:${lang}Date`,
        {'input_datetime': now.add(1, 'day')});
      expect(result, now);
    });

    test('Dataframe input/output', async () => {
      const df = DataFrame.fromCsv(`id,date,name\nid1,${Date.now()},datagrok`);
      const result = await grok.functions.call(`CVMTests:${lang}Dataframe`,
        {'df': df, 'dfNumerical': df, 'dfCategorical': df});
      expectTable(result['resultDf'], df);
      expectTable(result['resultNumerical'], df);
      expectTable(result['resultCategorical'], df);
    });

    test('Dataframe performance test 15 consequently', async () => {
      const results = [];
      for (let i = 0; i < 15; i++) {
        results.push(await getScriptTime(`CVMTests:${lang}SingleDf`,
          {'df': grok.data.demo.demog(10000)}));
      }
      const sum = results.reduce((p, c) => p + c, 0);
      return {'Average time': sum / results.length,
        'Min time': Math.min(...results), 'Max time': Math.max(...results)};
    }, {timeout: 120000});

    test('Dataframe performance test 25 parallel', async () => {
      const calls = [];
      for (let i = 0; i < 25; i++) {
        calls.push(getScriptTime(`CVMTests:${lang}SingleDf`,
          {'df': grok.data.demo.demog(10000)}));
      }
      const results = await Promise.all(calls);
      const sum = results.reduce((p, c) => p + c, 0);
      return {'Average time': sum / results.length,
        'Min time': Math.min(...results), 'Max time': Math.max(...results)};
    }, {timeout: 120000});

    test('Map type input/output', async () => {
      const result = await grok.functions.call(`CVMTests:${lang}Map`,
        {'input_map': {'hello': 'world'}, 'unique_key': 'my_key'});
      expectObject(result, {'hello': 'world', 'my_key': 'Datagrok'});
    });

    if (lang != 'Node' && lang !== 'JavaScript') {
      test('Graphics output, Column input', async () => {
        const df = DataFrame.fromCsv('x,y\n1,2\n3,4\n5,6');
        const result = await grok.functions.call(`CVMTests:${lang}Graphics`,
          {'df': df, 'xName': 'x', 'yName': 'y'});
        expect(result && result.length > 0, true);
      });
    }

    test('Column list', async () => {
      const df = DataFrame.fromCsv(`id,date,name\nid1,${Date.now()},datagrok`);
      const result = await grok.functions.call(`CVMTests:${lang}ColumnList`,
        {'df': df, 'cols': df.columns.toList()});
      df.columns.remove('id');
      expectTable(result, df);
    });

    test('Calculated column test', async () => {
      const df = DataFrame.fromCsv('x,y\n1,2\n3,4\n5,6');
      const column = await df.columns.addNewCalculated('new', `CVMTests:${lang}CalcColumn(\${x} + \${y})`);
      expect(df.columns.contains(column.name), true);
      expect(column.get(0), 6);
    });

    test('Calculated column performance', async () => {
      const df = grok.data.demo.demog(10000);
      const start = Date.now();
      await df.columns.addNewCalculated('new', `CVMTests:${lang}CalcColumn(\${age})`);
      return `Execution time: ${Date.now() - start}`;
    });

    test('File input', async () => {
      const files = await grok.dapi.files.list('System:AppData/CvmTests/', false, 'cars.csv');
      const result = await grok.functions.call(`CVMTests:${lang}LinesCount`,
        {'file': files[0], 'header': true, 'separator': ',', 'dec': '.'});
      const expected = ['Python', 'Octave'].includes(lang) ? 31 : 30;
      expect(result, expected);
    });

    if (lang === 'Python') {
      test('Environment string', async () => {
        const result = await grok.functions.call('CVMTests:PythonAnchorsCount',
          {'html': '  <ul>\n' +
                            '      <li><a href="https://example.com">Website</a></li>\n' +
                            '  <li><a href="mailto:m.bluth@example.com">Email</a></li>\n' +
                            '  <li><a href="tel:+123456789">Phone</a></li>\n' +
                            '  </ul>'});
        expect(result, 3);
      }, {timeout: 120000});

      test('File type input and environment yaml', async () => {
        const files = await grok.dapi.files.list('System:AppData/CvmTests/images', false, 'silver.jpg');
        const result = await grok.functions.call('CVMTests:ImagePixelCount',
          {'fileInput': files[0]});
        expect(79498, result);
      });
    }
  });
}

export async function getScriptTime(name: string, params: object = {}): Promise<number> {
  const start = Date.now();
  await grok.functions.call(name, params);
  return Date.now() - start;
}

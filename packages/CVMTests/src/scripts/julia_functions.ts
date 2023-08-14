import * as grok from 'datagrok-api/grok';
import {category, expect, expectObject, expectTable, test} from '@datagrok-libraries/utils/src/test';
import {DataFrame, FileSource} from 'datagrok-api/dg';
import dayjs from 'dayjs';

category('Scripts: Julia scripts', async () => {
  test('int, double, bool, string input/output', async () => {
    const int = 2;
    const double = 0.3;
    const bool = true;
    const str = 'Datagrok';
    const result = await grok.functions.call('CVMTests:JuliaSimple',
      {'integer_input': int, 'double_input': double, 'bool_input': bool, 'string_input': str});
    expectObject(result, {'integer_output': int, 'double_output': double,
      'bool_output': bool, 'string_output': str});
  }, {skipReason: 'Skipped for a while'});
});

import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import dayjs from 'dayjs';

const resCsv = `col1,new_col
val1,1
val2,1`;

category('Pyodide: FuncCall execution', async () => {
  test('Simple types call js', async () => {
    const func = DG.Func.byName('Pyodide:SimpleInputsCallJS');
    const fc = func.prepare({in1: true, in2: 2, in3: 3.5, in4: 'my string'});
    await fc.call();
    const {out1, out2, out3, out4} = fc.outputs;
    const actual = {out1, out2, out3, out4};
    const expected  = {out1: false, out2: 3, out3: 4.5, out4: 'my string edited'};
    expectDeepEqual(actual, expected);
  }, { skipReason: 'Old CI chrome version' });

  test('Simple types call pyodide', async () => {
    const func = DG.Func.byName('Pyodide:SimpleInputsCallPy');
    const fc = func.prepare({in1: true, in2: 2, in3: 3.5, in4: 'my string'});
    await fc.call();
    const {out1, out2, out3, out4} = fc.outputs;
    const actual = {out1, out2, out3, out4};
    const expected  = {out1: false, out2: 3, out3: 4.5, out4: 'my string edited'};
    expectDeepEqual(actual, expected);
  }, { skipReason: 'Old CI chrome version' });

  test('Datetime call js', async () => {
    const func = DG.Func.byName('Pyodide:DatetimeInputCallJS');
    const initial = dayjs('01-01-2020Z');
    const final = initial.add(1, 'day');
    const fc = func.prepare({input_datetime: initial});
    await fc.call();
    const {output_datetime} = fc.outputs;
    expectDeepEqual(output_datetime.unix(), final.unix());
  }, { skipReason: 'Old CI chrome version' });

  test('Datetime call pyodide', async () => {
    const func = DG.Func.byName('Pyodide:DatetimeInputCallPy');
    const initial = dayjs('01-01-2020Z');
    const final = initial.add(1, 'day');
    const fc = func.prepare({input_datetime: initial});
    await fc.call();
    const {output_datetime} = fc.outputs;
    expectDeepEqual(output_datetime.unix(), final.unix());
  }, { skipReason: 'Old CI chrome version' });

  test('Dataframe call js', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromList('string', 'col1', ['val1', 'val2'])]);
    const func = DG.Func.byName('Pyodide:DataFrameInputCallJS');
    const fc = func.prepare({input_df: df});
    await fc.call();
    const {output_df} = fc.outputs;
    expectDeepEqual(output_df, DG.DataFrame.fromCsv(resCsv));
  }, { skipReason: 'Old CI chrome version' });

  test('Dataframe call pyodide', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromList('string', 'col1', ['val1', 'val2'])]);
    const func = DG.Func.byName('Pyodide:DataFrameInputCallPy');
    const fc = func.prepare({input_df: df});
    await fc.call();
    const {output_df} = fc.outputs;
    expectDeepEqual(output_df, DG.DataFrame.fromCsv(resCsv));
  }, { skipReason: 'Old CI chrome version' });
});

import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {customDeepEqual} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';

function df(cols: Record<string, number[]>) {
  return DG.DataFrame.fromColumns(
    Object.entries(cols).map(([name, values]) => DG.Column.fromList('double', name, values)));
}

category('ComputeUtils: Driver customDeepEqual', async () => {
  test('Equal dataframes compare equal', async () => {
    expectDeepEqual(customDeepEqual(df({x: [1, 2], y: [3, 4]}), df({x: [1, 2], y: [3, 4]})), true);
  });

  test('Different values compare not equal', async () => {
    expectDeepEqual(customDeepEqual(df({x: [1, 2]}), df({x: [1, 5]})), false);
  });

  test('Same column count with different names compares not equal', async () => {
    expectDeepEqual(customDeepEqual(df({x: [1, 2]}), df({y: [1, 2]})), false);
  });
});

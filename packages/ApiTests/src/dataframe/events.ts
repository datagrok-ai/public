import * as DG from 'datagrok-api/dg';
import {category, expect, test} from "@datagrok-libraries/utils/src/test";

category('DataFrame', () => {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.TYPE.FLOAT, 'x', [1, 2, 3]),
    DG.Column.fromList(DG.TYPE.FLOAT, 'y', [4, 5, 6]),
    DG.Column.fromList(DG.TYPE.FLOAT, 'z', [7, 8, 9]),
  ]);
    
  test('raise an event on set', async () => {
    let dataChanged = false;
    df.onDataChanged.subscribe(() => {
      dataChanged = true;
    })
    df.set('x', 1, 0)
    expect(dataChanged, true);
  });

  test('raise an event on setValues', async () => {
    let dataChanged = false;
    df.onDataChanged.subscribe(() => {
      dataChanged = true;
    })
    df.rows.setValues(1, [0,0,0])
    expect(dataChanged, true);
  });
});
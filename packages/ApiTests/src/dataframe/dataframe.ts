import * as DG from 'datagrok-api/dg';
import {category, test} from "@datagrok-libraries/utils/src/test";

category('DataFrame', () => {
  test('create from arrays', async () => {
    let t = DG.DataFrame.fromColumns([
      DG.Column.fromList('int', 'int', [1, 2, 3]),
      DG.Column.fromList('double', 'double', [1.1, 2.1, 3.1]),
      DG.Column.fromList('string', 'string', ["a", "b", "c"]),
      DG.Column.fromList('object', 'object', [{}, null, {a: 1, b: 2}])
    ]);
  });

  test('create from columns', async () => {
    let t = DG.DataFrame.create(3);
    t.columns.add(DG.Column.fromStrings('countries', ['USA', 'Canada', 'Mexico']));
    t.columns.add(DG.Column.fromStrings('population', ['321', '35', '121']));
  });

  test('create from csv', async () => {
    let table = DG.DataFrame.fromCsv(
      `make, model,    cylinders, volume, price
Honda, Civic,    4,         1.4,    15000
Honda, Accord,   6,         1.8,    20000
BMW,   328i,     4,         1.7,    60000        
BMW,   535i,     6,         1.5,    35000
Tesla, Roadster, ,          1.6,    100000
Tesla, Model S,  ,          1.6,    120000`);
  });

  test('create from json', async () => {
    let table = DG.DataFrame.fromJson(`[
  {
    "name": "Roger Federer",
    "height": 185,
    "born": "August 8, 1981"
  },
  {
    "name": "Rafael Nadal",
    "height": 185,
    "born": "June 3, 1986"
  }
]`);
  });

  test('create from typed arrays', async () => {
    let table = DG.DataFrame.fromColumns([
      DG.Column.fromInt32Array('ints', new Int32Array(100)),
      DG.Column.fromFloat32Array('floats', new Float32Array(100))
    ]);
  });

});
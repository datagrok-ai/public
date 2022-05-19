import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, expect, test} from "@datagrok-libraries/utils/src/test";

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

  //creation of dataframes used in testing
  let df1 = DG.DataFrame.create(2);
  df1.columns.add(DG.Column.fromStrings('countries', ['USA', 'Canada']));
  df1.columns.add(DG.Column.fromInt32Array('population', Int32Array.from([1, 4])));
  let df2 = DG.DataFrame.create(2);
  df2.columns.add(DG.Column.fromStrings('countries', ['France', 'Mexico']));
  df2.columns.add(DG.Column.fromInt32Array('population', Int32Array.from([2, 3])));
  let df = DG.DataFrame.create(4);
  df.columns.add(DG.Column.fromStrings('countries', ['USA', 'Canada', 'France', 'Mexico']));
  df.columns.add(DG.Column.fromInt32Array('population', Int32Array.from([1, 4, 2, 3])));
  
  test('append method', async () => {
    df1.append(df2);
    //expect(df1.append(df2), df);
  });

  test('cell method', async () => {
    //df.cell(1, 'countries');
    expect(df1.cell(1, 'countries').toString(), 'countries : 1');
  });

  test('get method', async () => {
    expect(df.get('countries', 1), 'Canada');
  });

  test('change column type method', async () => {
    df1.changeColumnType('countries', 'int');
    expect(typeof(df1.get('countries', 1)), 'number');
  });

  test('col method', async () => {
    let res = df1.col('population')?.getRawData();
    expect(res?.toString(), "1,4");
  });

  test('getcol method', async () => {
    let res = df1.getCol('population').getRawData();
    expect(res.toString(), "1,4");
  });

  test('get sorted order method', async () => {
    const arr: object[] = [];
    arr.push(df.columns.names()[1]);
    let order = df.getSortedOrder(arr);
    expect(order.toString(), '0,2,3,1');
  });

  test('get table info', async () => {
    return df1.getTableInfo();
  });

  test('toString method', async () => {
    expect(typeof(df1.getCol('population').getRawData().toString()), 'string');
  });

  test('toCsv method', async () => {
    return df1.toCsv();
  });

  test('set method', async () => {
    df1.set("population", 1, 5);
    expect(5, df1.get("population", 1));
  });

  test('get density method', async () => {
    let t = DG.DataFrame.create(3);
    t.columns.add(DG.Column.fromInt32Array('numbers', Int32Array.from([1, 4, 2])));
    t.columns.add(DG.Column.fromInt32Array('digits', Int32Array.from([2, 1, 3])));
    let res = t.getDensity(4, 2, 'numbers', 'digits');
    return res;
  });

  test('set tag method and get tag method', async () => {
    df1.columns.byName('population').setTag('units', 'm');
    expect('m', df1.columns.byName('population').getTag('units'));
  });

  //test('row method', async () => {
    //return df1.row(1);
  //});

});
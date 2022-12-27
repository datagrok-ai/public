import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import dayjs from 'dayjs';
import {hashDataFrame} from '@datagrok-libraries/utils/src/dataframe-utils';

category('DataFrame', () => {
  const createDf = (): DG.DataFrame => {
    const df = DG.DataFrame.create(2);
    df.columns.add(DG.Column.fromStrings('countries', ['USA', 'Canada']));
    df.columns.add(DG.Column.fromInt32Array('population', Int32Array.from([1, 4])));
    return df;
  };

  const createDf2 = (): DG.DataFrame => {
    const df = DG.DataFrame.create(4);
    df.columns.add(DG.Column.fromStrings('countries', ['USA', 'Canada', 'France', 'Mexico']));
    df.columns.add(DG.Column.fromInt32Array('population', Int32Array.from([1, 4, 2, 3])));
    return df;
  };

  test('byte array', async () => {
    const t = grok.data.testData('demog');
    //const data = t.toByteArray();
    //const t2 = DG.DataFrame.fromByteArray(data);
    //expect(t.toCsv(), t2.toCsv());
  });

  test('create from arrays', async () => {
    DG.DataFrame.fromColumns([
      DG.Column.fromList('int', 'int', [1, 2, 3]),
      DG.Column.fromList('double', 'double', [1.1, 2.1, 3.1]),
      DG.Column.fromList('string', 'string', ['a', 'b', 'c']),
      DG.Column.fromList('object', 'object', [{}, null, {a: 1, b: 2}]),
    ]);
  });

  test('create from columns', async () => {
    const t = DG.DataFrame.create(3);
    t.columns.add(DG.Column.fromStrings('countries', ['USA', 'Canada', 'Mexico']));
    t.columns.add(DG.Column.fromStrings('population', ['321', '35', '121']));
  });

  test('create from csv', async () => {
    DG.DataFrame.fromCsv(
      `make, model,    cylinders, volume, price
Honda, Civic,    4,         1.4,    15000
Honda, Accord,   6,         1.8,    20000
BMW,   328i,     4,         1.7,    60000        
BMW,   535i,     6,         1.5,    35000
Tesla, Roadster, ,          1.6,    100000
Tesla, Model S,  ,          1.6,    120000`);
  });

  test('create from json', async () => {
    DG.DataFrame.fromJson(`[
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
    DG.DataFrame.fromColumns([
      DG.Column.fromInt32Array('ints', new Int32Array(100)),
      DG.Column.fromFloat32Array('floats', new Float32Array(100)),
    ]);
  });

  test('append method', async () => {
    const df1 = createDf();
    const df2 = createDf();
    df1.append(df2);
    //expect(df1.append(df2), df);
  });

  test('cell method', async () => {
    //df.cell(1, 'countries');
    const df1 = createDf();
    expect(df1.cell(1, 'countries').toString(), 'countries : 1');
  });

  test('get method', async () => {
    const df = createDf2();
    expect(df.get('countries', 1), 'Canada');
  });

  test('change column type method', async () => {
    const df1 = createDf();
    df1.changeColumnType('countries', 'int');
    expect(typeof(df1.get('countries', 1)), 'number');
  });

  test('col method', async () => {
    const df1 = createDf();
    const res = df1.col('population')?.getRawData();
    expect(res?.toString(), '1,4');
  });

  test('getcol method', async () => {
    const df1 = createDf();
    const res = df1.getCol('population').getRawData();
    expect(res.toString(), '1,4');
  });

  test('get sorted order method', async () => {
    const df = createDf2();
    const arr = [];
    arr.push(df.columns.names()[1]);
    const order = df.getSortedOrder(arr);
    expect(order.toString(), '0,2,3,1');
  });

  test('get table info', async () => {
    const df1 = createDf();
    return df1.getTableInfo();
  });

  test('toString method', async () => {
    const df1 = createDf();
    expect(typeof(df1.getCol('population').getRawData().toString()), 'string');
  });

  test('toCsv method', async () => {
    const df1 = createDf();
    return df1.toCsv();
  });

  test('set method', async () => {
    const df1 = createDf();
    df1.set('population', 1, 5);
    expect(5, df1.get('population', 1));
  });

  test('DataFrame.getDensity', async () => {
    const t = DG.DataFrame.create(3);
    t.columns.add(DG.Column.fromInt32Array('numbers', Int32Array.from([1, 4, 2])));
    t.columns.add(DG.Column.fromInt32Array('digits', Int32Array.from([2, 1, 3])));
    return t.getDensity(4, 2, 'numbers', 'digits');
  });

  test('set tag method and get tag method', async () => {
    const df1 = createDf();
    df1.columns.byName('population').setTag('units', 'm');
    expect('m', df1.columns.byName('population').getTag('units'));
  });

  test('column list add', async () => {
    const df1 = createDf();
    df1.columns.add(DG.Column.fromInt32Array('popularity', Int32Array.from([12, 10])));
    const res = df1.columns.byName('popularity').getRawData();
    expect(res.toString(), '12,10');
  });

  test('column list addNew', async () => {
    const df = createDf2();
    df.columns.addNew('newColumn', 'string');
    expect(typeof(df.get('newColumn', 1)), 'string');
  });

  test('column list addNewBool', async () => {
    const df = createDf2();
    df.columns.addNewBool('newColumnBool');
    expect(typeof(df.get('newColumnBool', 1)), 'boolean');
  });

  test('column list addNewDateTime', async () => {
    const df = createDf2();
    df.columns.addNewDateTime('newColumnDateTime');
    expect(typeof(df.get('newColumnDateTime', 1)), 'object');
  });

  test('column list addNewFloat', async () => {
    const df = createDf2();
    df.columns.addNewFloat('newColumnFloat');
    expect(typeof(df.get('newColumnFloat', 1)), 'number');
  });

  test('column list addNewInt', async () => {
    const df = createDf2();
    df.columns.addNewInt('newColumnInt');
    expect(typeof(df.get('newColumnInt', 1)), 'number');
  });

  test('column list addNewQnum', async () => {
    const df = createDf2();
    df.columns.addNewQnum('newColumnQnum');
    expect(typeof(df.get('newColumnQnum', 1)), 'number');
  });

  test('column list addNewString', async () => {
    const df = createDf2();
    df.columns.addNewString('newColumnString');
    expect(typeof(df.get('newColumnString', 1)), 'string');
  });

  test('column list byIndex', async () => {
    const df = createDf2();
    expect(df.columns.byIndex(1).getRawData().toString(), '1,4,2,3');
  });

  test('column list byName', async () => {
    const df = createDf2();
    expect(df.columns.byName('population').getRawData().toString(), '1,4,2,3');
  });

  test('column list bySemType', async () => {
    const df1 = createDf();
    df1.onSemanticTypeDetected.subscribe((_) => {
      const res = df1.columns.bySemType('country');
      expect(res?.toString(), 'countries');
    });
  });

  test('column list bySemTypeAll', async () => {
    const df1 = createDf();
    df1.onSemanticTypeDetected.subscribe((_) => {
      const res = df1.columns.bySemTypeAll('country');
      expect(res.toString(), 'countries');
    });
  });

  test('column list contains', async () => {
    const df = createDf2();
    expect(df.columns.contains('population'), true);
  });

  test('column list getUnusedName', async () => {
    const df = createDf2();
    expect(df.columns.getUnusedName('population').toString(), 'population (2)');
  });

  test('column list insert', async () => {
    const df1 = createDf();
    const newColumn = DG.Column.fromStrings('data', ['12', '34']);
    df1.columns.insert(newColumn);
    expect(df1.get('data', 1).toString(), '34');
  });

  test('column list names', async () => {
    const df3 = createDf();
    expect(df3.columns.names().toString(), 'countries,population');
  });

  test('column list remove', async () => {
    const df3 = createDf();
    df3.columns.remove('population');
    expect(df3.columns.names().toString(), 'countries');
  });

  test('column list replace', async () => {
    const df4= createDf();
    const newColumn = DG.Column.fromStrings('data', ['12', '34']);
    df4.columns.replace('countries', newColumn);
    expect(df4.columns.names().toString(), 'data,population');
  });

  test('column list toList', async () => {
    const df1 = createDf();
    return df1.columns.toList();
  });

  test('column list toString', async () => {
    const df1 = createDf();
    return df1.columns.toString();
  });

  test('row list addNew', async () => {
    const df5 = createDf();
    df5.rows.addNew(['12', 23]);
    expect(df5.get('countries', df5.rowCount - 1).toString(), '12');
  });

  test('row list filter', async () => {
    const df1 = createDf();
    return df1.rows.filter((row) => row.countries === 'USA');
  });

  test('row list insertAt', async () => {
    const df6 = createDf();
    df6.rows.insertAt(2, 2);
    expect(df6.get('countries', 2), '');
  });

  test('row list match', async () => {
    const df5 = createDf();
    return df5.rows.match('countries = USA').highlight();
  });

  test('datetime column', async () => {
    const t = grok.data.testData('demog');
    const c = t.columns.byName('started');
    c.set(1, dayjs('2022-01-01'));
    expect(c.get(1).valueOf(), 1640984400000);
    c.set(1, null);
    expect(c.get(1), null);
    const v = grok.shell.addTableView(t);
    v.close();
  });

  test('hash', async () => {
    const df = grok.data.demo.demog(100000);
    expect(hashDataFrame(df).length, 32);

    const df1 = DG.DataFrame.fromCsv(`a,b\n1,0\n2,0\n3,0`);
    const df2 = DG.DataFrame.fromCsv(`a,b\n2,0\n1,0\n3,0`);
    expectArray(hashDataFrame(df1), hashDataFrame(df2));

    const df3 = DG.DataFrame.fromCsv(`a,b\n"abc",0\n"dce",0\n"xyz",0`);
    const df4 = DG.DataFrame.fromCsv(`a,b\n"dce",0\n"abc",0\n"xyz",0`);
    expectArray(hashDataFrame(df3), hashDataFrame(df4));
  });
});

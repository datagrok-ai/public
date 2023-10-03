import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import dayjs from 'dayjs';
import {hashDataFrame} from '@datagrok-libraries/utils/src/dataframe-utils';
import {category, expect, expectArray, test, expectTable} from '@datagrok-libraries/utils/src/test';

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

category('DataFrame: Methods', () => {

  test('create from typed arrays', async () => {
    DG.DataFrame.fromColumns([
      DG.Column.fromInt32Array('ints', new Int32Array(100)),
      DG.Column.fromFloat32Array('floats', new Float32Array(100)),
    ]);
  });

  test('append', async () => {
    const df1 = createDf();
    const df2 = DG.DataFrame.create(2);
    df2.columns.add(DG.Column.fromStrings('countries', ['France', 'Mexico']));
    df2.columns.add(DG.Column.fromInt32Array('population', Int32Array.from([2, 3])));
    df1.append(df2, true);
    expectTable(df1, createDf2());
  });

  test('cell', async () => {
    const df1 = createDf();
    expect(df1.cell(1, 'countries').toString(), 'countries : 1');
  });

  test('changeColumnType', async () => {
    const df1 = createDf();
    df1.changeColumnType('countries', 'int');
    expect(typeof(df1.get('countries', 1)), 'number');
  });

  test('clone', async () => {
    const df = createDf2();
    const dfCloned = df.clone();
    expectTable(dfCloned, df);
  });

  test('col', async () => {
    const df1 = createDf();
    const res = df1.col('population')?.getRawData();
    expect(res?.toString(), '1,4');
  });

  test('fireValuesChanged', async () => {
    const df = createDf();
    let fired = false;
    const sub = df.onValuesChanged.subscribe(() => fired = true);
    df.fireValuesChanged();
    sub.unsubscribe();
    expect(fired, true);
  });

  test('get', async () => {
    const df = createDf2();
    expect(df.get('countries', 1), 'Canada');
  });

  test('getCol', async () => {
    const df1 = createDf();
    const res = df1.getCol('population').getRawData();
    expect(res.toString(), '1,4');
  });

  test('getDensity', async () => {
    const t = DG.DataFrame.create(3);
    t.columns.add(DG.Column.fromInt32Array('numbers', Int32Array.from([1, 4, 2])));
    t.columns.add(DG.Column.fromInt32Array('digits', Int32Array.from([2, 1, 3])));
    return t.getDensity(4, 2, 'numbers', 'digits');
  });

  test('getSortedOrder', async () => {
    const df = createDf2();
    const arr = [];
    arr.push(df.columns.names()[1]);
    const order = df.getSortedOrder(arr);
    expect(order.toString(), '0,2,3,1');
  });

  test('getTableInfo', async () => {
    const df1 = createDf();
    return df1.getTableInfo();
  });

  test('getTag | setTag', async () => {
    const df1 = createDf();
    df1.columns.byName('population').setTag('units', 'm');
    expect('m', df1.columns.byName('population').getTag('units'));
  });

  test('groupBy', async () => {
    const df = createDf();
    const df1 = DG.DataFrame.create(2);
    df1.columns.add(DG.Column.fromStrings('countries', ['Canada', 'USA']));
    expectTable(df.groupBy(['countries']).aggregate(), df1);
  });

  test('row', async () => {
    const df = createDf();
    expectArray([...df.row(0).cells].map((i) => i.value), ['USA', 1]);
  });

  test('set', async () => {
    const df1 = createDf();
    df1.set('population', 1, 5);
    expect(5, df1.get('population', 1));
  });

  test('toByteArray | fromByteArray', async () => {
    const t = createDf();
    const data = t.toByteArray();
    const t2 = DG.DataFrame.fromByteArray(data);
    expectTable(t2, t);
  });

  test('toCsv | fromCsv', async () => {
    const t = createDf();
    const data = t.toCsv();
    const t2 = DG.DataFrame.fromCsv(data);
    expectTable(t2, t);
  });

  test('toJson | fromJson', async () => {
    const t = createDf();
    //@ts-ignore
    const data = JSON.stringify(t.toJson());
    const t2 = DG.DataFrame.fromJson(data);
    expectTable(t2, t);
  });
  
  test('toString', async () => {
    const df = createDf();
    df.name = 'Table';
    expect(df.toString(), 'Table (2 rows, 2 columns)');
  });

  test('unpivot', async () => {
    const df = createDf();
    const df1 = DG.DataFrame.create(2);
    df1.columns.add(DG.Column.fromStrings('countries', ['USA', 'Canada']));
    df1.columns.add(DG.Column.fromStrings('Category', ['population', 'population']));
    df1.columns.add(DG.Column.fromInt32Array('Value', Int32Array.from([1, 4])));
    expectTable(df.unpivot(['countries'], ['population']), df1);
  });

  test('create', async () => {
    const df = createDf();
    const str = `countries,population
USA,1
Canada,4`;
    expect(df.toCsv(), str);
  });

  test('fromColumns', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList('int', 'int', [1, 2, 3]),
      DG.Column.fromList('double', 'double', [1.1, 2.1, 3.1]),
      DG.Column.fromList('string', 'string', ['a', 'b', 'c']),
      DG.Column.fromList('object', 'object', [{}, null, {a: 1, b: 2}]),
    ]);
    const str = `int,double,string,object
1,1.100000023841858,a,[object Object]
2,2.0999999046325684,b,
3,3.0999999046325684,c,[object Object]`;
    expect(df.toCsv(), str);
  });

  test('fromObjects', async () => {
    const t = createDf();
    //@ts-ignore
    const data = t.toJson();
    const t2 = DG.DataFrame.fromObjects(data);
    expectTable(t2!, t);
  });

  test('fromProperties', async () => {
    
  });
});

category('DataFrame: ColumnList', () => {
  test('create from columns', async () => {
    const t = DG.DataFrame.create(3);
    t.columns.add(DG.Column.fromStrings('countries', ['USA', 'Canada', 'Mexico']));
    t.columns.add(DG.Column.fromStrings('population', ['321', '35', '121']));
  });

  test('add', async () => {
    const df1 = createDf();
    df1.columns.add(DG.Column.fromInt32Array('popularity', Int32Array.from([12, 10])));
    const res = df1.columns.byName('popularity').getRawData();
    expect(res.toString(), '12,10');
  });

  test('addNew', async () => {
    const df = createDf2();
    df.columns.addNew('newColumn', 'string');
    expect(typeof(df.get('newColumn', 1)), 'string');
  });

  test('addNewBool', async () => {
    const df = createDf2();
    df.columns.addNewBool('newColumnBool');
    expect(typeof(df.get('newColumnBool', 1)), 'boolean');
  });

  test('addNewDateTime', async () => {
    const df = createDf2();
    df.columns.addNewDateTime('newColumnDateTime');
    expect(typeof(df.get('newColumnDateTime', 1)), 'object');
  });

  test('addNewFloat', async () => {
    const df = createDf2();
    df.columns.addNewFloat('newColumnFloat');
    expect(typeof(df.get('newColumnFloat', 1)), 'number');
  });

  test('addNewInt', async () => {
    const df = createDf2();
    df.columns.addNewInt('newColumnInt');
    expect(typeof(df.get('newColumnInt', 1)), 'number');
  });

  test('addNewQnum', async () => {
    const df = createDf2();
    df.columns.addNewQnum('newColumnQnum');
    expect(typeof(df.get('newColumnQnum', 1)), 'number');
  });

  test('addNewString', async () => {
    const df = createDf2();
    df.columns.addNewString('newColumnString');
    expect(typeof(df.get('newColumnString', 1)), 'string');
  });

  test('byIndex', async () => {
    const df = createDf2();
    expect(df.columns.byIndex(1).getRawData().toString(), '1,4,2,3');
  });

  test('byName', async () => {
    const df = createDf2();
    expect(df.columns.byName('population').getRawData().toString(), '1,4,2,3');
  });

  test('bySemType', async () => {
    const df1 = createDf();
    df1.onSemanticTypeDetected.subscribe((_) => {
      const res = df1.columns.bySemType('country');
      expect(res?.toString(), 'countries');
    });
  });

  test('bySemTypeAll', async () => {
    const df1 = createDf();
    df1.onSemanticTypeDetected.subscribe((_) => {
      const res = df1.columns.bySemTypeAll('country');
      expect(res.toString(), 'countries');
    });
  });

  test('contains', async () => {
    const df = createDf2();
    expect(df.columns.contains('population'), true);
  });

  test('getUnusedName', async () => {
    const df = createDf2();
    expect(df.columns.getUnusedName('population').toString(), 'population (2)');
  });

  test('insert', async () => {
    const df1 = createDf();
    const newColumn = DG.Column.fromStrings('data', ['12', '34']);
    df1.columns.insert(newColumn);
    expect(df1.get('data', 1).toString(), '34');
  });

  test('names', async () => {
    const df3 = createDf();
    expect(df3.columns.names().toString(), 'countries,population');
  });

  test('remove', async () => {
    const df3 = createDf();
    df3.columns.remove('population');
    expect(df3.columns.names().toString(), 'countries');
  });

  test('replace', async () => {
    const df4= createDf();
    const newColumn = DG.Column.fromStrings('data', ['12', '34']);
    df4.columns.replace('countries', newColumn);
    expect(df4.columns.names().toString(), 'data,population');
  });

  test('toList', async () => {
    const df1 = createDf();
    expect(df1.columns.toList().length, 2);
  });

  test('toString', async () => {
    const df1 = createDf();
    return df1.columns.toString();
  });
});

category('DataFrame: RowList', () => {
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

  test('row list removeWhere', async () => {
    const df = createDf2();
    df.rows.removeWhere((r) => {
      return r.get('population') < 3;
    });
    expect(df.rowCount, 2);
  });

  test('datetime column', async () => {
    const t = grok.data.testData('demog');
    const c = t.columns.byName('started');
    c.set(1, dayjs.utc('2022-01-01'));
    expect(c.get(1).valueOf(), 1640995200000);
    c.set(1, null);
    expect(c.get(1), null);
  });

  test('hash', async () => {
    const df = grok.data.demo.demog(10);
    expect(hashDataFrame(df).length, 32);

    const df1 = DG.DataFrame.fromCsv(`a,b\n1,0\n2,0\n3,0`);
    const df2 = DG.DataFrame.fromCsv(`a,b\n2,0\n1,0\n3,0`);
    expectArray(hashDataFrame(df1), hashDataFrame(df2));

    const df3 = DG.DataFrame.fromCsv(`a,b\n"abc",0\n"dce",0\n"xyz",0`);
    const df4 = DG.DataFrame.fromCsv(`a,b\n"dce",0\n"abc",0\n"xyz",0`);
    expectArray(hashDataFrame(df3), hashDataFrame(df4));
  });

  test('emptyDataFrameToCsv', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromColumns([]);
    const csv: string = df.toCsv();
    expect(csv, '');
  });
});

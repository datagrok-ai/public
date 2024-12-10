import * as DG from 'datagrok-api/dg';
// import * as grok from 'datagrok-api/grok';

import {category, expect, expectArray, test, expectTable, expectFloat} from '@datagrok-libraries/utils/src/test';

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

const createCol = (): DG.Column => DG.Column.fromStrings('countries', ['USA', 'Canada', 'France', 'Mexico']);

const createCol2 = (): DG.Column => DG.Column.fromInt32Array('population', Int32Array.from([1, 4, 2, 3]));

// DataFrame

category('DataFrame: Methods', () => {

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
    t.getDensity(4, 2, 'numbers', 'digits');
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
    df1.getTableInfo();
  });

  test('getTag | setTag', async () => {
    const df1 = createDf();
    df1.columns.byName('population').meta.units = 'm';
    expect('m', df1.columns.byName('population').meta.units);
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
    expect(t2 instanceof DG.DataFrame);
    expectTable(t2, t);
    const df = DG.DataFrame.fromColumns([]);
    const csv = df.toCsv();
    expect(csv, '');
  });

  test('toCsv with grid settings', async () => {
    const t = DG.DataFrame.fromCsv(`x, y, z, a
      1, 6, 3, 4
      5, 2, 7, 8
      9, 10, 11, 12`);
    const grid = DG.Grid.create(t);
    grid.columns.byName('a')!.visible = false;
    grid.columns.setOrder(['x', 'z', 'y', 'a']);
    grid.sort(['x'], [false]);

    const t2 = DG.DataFrame.fromCsv(t.toCsv({visibleColumnsOnly: true}, grid));
    expect(t2.columns.length === 3);
    expect(t2.columns.names()[1] === 'z');
    expect(t2.columns.names()[2] === 'y');
    expect(t2.columns.byName('x').getNumber(0) === 9);
    expect(t2.columns.byName('z').getNumber(1) === 7);
    expect(t2.columns.byName('y').getNumber(2) === 6);
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
    const prop1 = DG.Property.fromOptions({name: 'countries', type: 'string'});
    const prop2 = DG.Property.fromOptions({name: 'population', type: 'int'});
    const df = DG.DataFrame.fromProperties([prop1, prop2], 4);
    expect(df.columns.length, 2);
    expect(df.rowCount, 4);
  });
});

// Column

category('DataFrame: Column', () => {
  const COL1 = createCol();
  const COL2 = createCol2();

  test('aggregate', async () => {
    expect(COL1.aggregate(DG.STR_AGG.CONCAT_ALL), 'USA,Canada,France,Mexico');
    expect(COL2.aggregate(DG.AGG.SUM), 10);
  });

  test('applyFormula', async () => {
    const df = createDf2();
    const col = df.getCol('population');
    await col.applyFormula('${population} + 5', 'int');
    expectArray(df.getCol('population').toList(), [6, 9, 7, 8]);
  });

  // test('compact', async () => {
  // });

  test('convertTo', async () => {
    const col = createCol2();
    expectArray(col.convertTo('string').toList(), ['1', '4', '2', '3']);
  });

  test('get', async () => {
    expect(COL2.get(0), 1);
    expect(COL2.get(1), 4);
  });

  test('getCategory', async () => {
    expect(COL1.getCategory(0), 'Canada');
    expect(COL1.getCategory(1), 'France');
  });

  test('getCategoryOrder | setCategoryOrder', async () => {
    const col = createCol();
    col.setCategoryOrder(['Canada', 'France', 'USA', 'Mexico']);
    expectArray(col.getCategoryOrder(), ['Canada', 'France', 'USA', 'Mexico']);
  });

  test('getRawData | setRawData', async () => {
    const col = DG.Column.int('col', 3);
    col.setRawData(Int32Array.from([1, 2, 3]));
    expectArray(col.getRawData(), [1, 2, 3]);
  }, {skipReason:'GROK-16406'});

  test('getSortedOrder', async () => {
    expectArray(COL1.getSortedOrder(), [1, 2, 3, 0]);
    expectArray(COL2.getSortedOrder(), [0, 2, 3, 1]);
  });

  test('getString | setString', async () => {
    const col = createCol();
    expect(col.getString(0), 'USA');
    col.setString(1, '20');
    expect(col.getString(1), '20');
  });

  test('getTag | setTag', async () => {
    const col = createCol();
    col.setTag('tag', 'value');
    expect(col.getTag('tag'), 'value');
  });

  test('init', async () => {
    const col = DG.Column.int('col', 5);
    col.init((i) => i + 1);
    expectArray(col.toList(), [1, 2, 3, 4, 5]);
  });

  test('isNone', async () => {
    const col = DG.Column.int('col', 5);
    expect(col.isNone(0), true);
    expect(COL1.isNone(0), false);
  });

  test('matches', async () => {
    expect(COL1.matches('string'), true);
    expect(COL1.matches('int'), false);
  });

  test('scale', async () => {
    expect(COL2.scale(0), 0);
    expect(COL2.scale(1), 1);
    expectFloat(COL2.scale(2), 0.33, 0.01);
    expectFloat(COL2.scale(3), 0.66, 0.01);
  });

  test('set', async () => {
    const col = createCol();
    col.set(0, 'UK');
    expect(col.get(0), 'UK');
  });
  
  test('toList', async () => {
    expectArray(createCol2().toList(), [1, 4, 2, 3]);
  });
  
  test('toString', async () => {
    expect(COL1.toString(), 'countries');
  });
  
  test('values', async () => {
    const gen = COL1.values()[Symbol.iterator]();
    expect(gen.next().value, 'USA');
    expect(gen.next().value, 'Canada');
  });
  
  test('bool', async () => {
    const col = DG.Column.bool('col', 3);
    col.init((_) => true);
    expect(col.type, 'bool');
    expectArray(col.toList(), [true, true, true]);
  });
  
  test('dataFrame', async () => {
    const col = DG.Column.dataFrame('col', 3);
    expect(col.type, 'dataframe');
  });
  
  test('dateTime', async () => {
    const col = DG.Column.dateTime('col', 3);
    col.init((_) => Date.now());
    expect(col.type, 'datetime');
  });
  
  test('float', async () => {
    const col = DG.Column.float('col', 3);
    col.init((_) => 1.1);
    expect(col.type, 'double');
    expectFloat(col.get(0) ?? 0, 1.1);
  });
  
  // test('fromBitSet', async () => {
    
  // });
  
  test('fromFloat32Array', async () => {
    const col = DG.Column.fromFloat32Array('col', Float32Array.from([1.23, 4.45, 3.34]), 3);
    expectFloat(col.get(0) ?? 0, 1.23);
  });
  
  test('fromInt32Array', async () => {
    const col = DG.Column.fromInt32Array('col', Int32Array.from([1, 2, 3]), 3);
    expect(col.get(0) ?? 0, 1);
  });

  test('fromList', async () => {
    const col = DG.Column.fromList('string', 'col', ['a', 'b', 'c']);
    expect(col.get(0), 'a');
  });

  test('fromStrings', async () => {
    const col = DG.Column.fromStrings('col', ['a', 'b', 'c']);
    expect(col.get(0), 'a');
  });

  test('fromType', async () => {
    const col = DG.Column.fromType('bigint', 'col', 5);
    expect(col.type, 'bigint');
  });

  test('int', async () => {
    const col = DG.Column.int('col', 3);
    col.init((_) => 1);
    expect(col.get(1), 1);
  });

  test('qnum', async () => {
    const col = DG.Column.qnum('col', 3, [1, 2, 3]);
    expectFloat(col.get(0) ?? 0, 1);
  });

  test('string', async () => {
    const col = DG.Column.string('col', 3);
    col.init((_) => 'val');
    expect(col.get(0), 'val');
  });
});

// ColumnList

category('DataFrame: ColumnList', () => {

  test('iterator', async () => {
    const df = createDf();
    const gen = df.columns[Symbol.iterator]();
    const col: DG.Column = gen.next().value;
    expect(col.get(0), 'USA');
    expect(col.get(1), 'Canada');
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

  test('addNewBytes', async () => {
    const df = createDf();
    df.columns.addNewBytes('newColumnBytes');
    expect(typeof(df.get('newColumnBytes', 1)), 'object');
  });

  // addNewCalculated in separate file

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

  test('addNewVirtual', async () => {
    const df = createDf();
    df.columns.addNewVirtual('newColumnVirtual', (i) => i * 2, DG.TYPE.INT);
    expect(typeof(df.get('newColumnVirtual', 1)), 'number');
  });

  test('byIndex', async () => {
    const df = createDf2();
    expect(df.columns.byIndex(1).getRawData().toString(), '1,4,2,3');
  });

  test('byName', async () => {
    const df = createDf2();
    expect(df.columns.byName('population').getRawData().toString(), '1,4,2,3');
  });

  test('byNames', async () => {
    const df = createDf();
    const [col1, col2] = df.columns.byNames(['countries', 'population']);
    expect(col1.toList().toString(), 'USA,Canada');
    expect(col2.getRawData().toString(), '1,4');
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

  test('bySemTypesExact', async () => {
    const df1 = createDf();
    df1.onSemanticTypeDetected.subscribe((_) => {
      const res = df1.columns.bySemTypesExact(['country']);
      expect(res?.toString(), 'countries');
    });
  });

  test('byTags', async () => {
    const df = createDf();
    df.col('countries')?.setTag('tag1', 'value1');
    df.col('population')?.setTag('tag2', 'value2');
    const res1 = [...df.columns.byTags({tag1: 'value1'})];
    const res2 = [...df.columns.byTags({tag2: undefined})];
    expect(res1[0].name, 'countries');
    expect(res2[0].name, 'population');
  });

  test('contains', async () => {
    const df = createDf2();
    expect(df.columns.contains('population'), true);
  });

  test('getOrCreate', async () => {
    const df = createDf();
    df.columns.getOrCreate('countries', 'string');
    expect(df.columns.length, 2);
    df.columns.getOrCreate('countries1', 'string');
    expect(df.columns.length, 3);
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
    df1.columns.toString();
  });
});

// Row

category('DataFrame: Row', () => {
  const row = createDf().row(0);

  test('get', async () => {
    expect(row.get('countries'), 'USA');
    expect(row.get('population'), 1);
  });

  test('toDart', async () => {
    expect(typeof row.toDart(), 'object');
  });
});

// RowList

category('DataFrame: RowList', () => {

  test('iterator', async () => {
    const df = createDf();
    const gen = df.rows[Symbol.iterator]();
    const row = gen.next().value as DG.Row;
    expect(row.countries, 'USA');
    expect(row.population, 1);
  });

  test('addFilterState', async () => {
    const state = {
      column: 'population',
      type: 'string',
    };
    const df = createDf();
    df.rows.addFilterState(state);
  });

  test('addNew', async () => {
    const df5 = createDf();
    df5.rows.addNew(['12', 23]);
    expect(df5.get('countries', df5.rowCount - 1).toString(), '12');
  });

  test('filter', async () => {
    const df = createDf2();
    df.rows.filter((r) => r.get('population') < 3);
    expect(df.filter.trueCount, 2);
  });

  test('get', async () => {
    const df = createDf();
    expectArray([...df.rows.get(0).cells].map((i) => i.value), ['USA', 1]);
  });

  test('highlight', async () => {
    const df = createDf2();
    df.rows.highlight((i) => i > 1);
  });

  test('insertAt', async () => {
    const df6 = createDf();
    df6.rows.insertAt(2, 2);
    expect(df6.get('countries', 2), '');
  });

  test('match', async () => {
    const df5 = createDf();
    return df5.rows.match('countries = USA').highlight();
  });

  test('removeAt', async () => {
    const df = createDf2();
    df.rows.removeAt(1, 2);
    expect(df.rowCount, 2);
  });

  test('removeWhere', async () => {
    const df = createDf2();
    df.rows.removeWhere((r) => r.get('population') < 3);
    expect(df.rowCount, 2);
  });

  test('removeWhereIdx', async () => {
    const df = createDf2();
    df.rows.removeWhereIdx((i) => i > 1);
    expect(df.rowCount, 2);
  });

  test('requestFilter', async () => {
    const df = createDf2();
    df.rows.filter((r) => r.get('population') < 3);
    expect(df.filter.trueCount, 2);
    df.rows.requestFilter();
    expect(df.filter.trueCount, 4);
  });

  test('select', async () => {
    const df = createDf2();
    df.rows.select((r) => r.get('population') < 3);
    expect(df.selection.trueCount, 2);
  });

  test('setValues', async () => {
    const df = createDf();
    df.rows.setValues(0, ['France', 2]);
    expectArray([...df.rows.get(0).cells].map((i) => i.value), ['France', 2]);
  });

  test('toString', async () => {
    const str: string = createDf().rows.toString();
    expect(str.startsWith('(Instance of '), true);
  });
});

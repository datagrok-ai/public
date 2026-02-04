import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import { before, category, expect, test, expectArray, after } from '@datagrok-libraries/test/src/test';


category('Dapi: connection', () => {
  const dcParams = {
    dataSource: 'PostgresDart', server: 'localhost:5432', db: 'datagrok_dev', login: 'datagrok_dev', password: '123'
  };

  test('Create, save, delete, share', async () => {
    let dc = DG.DataConnection.create('Local DG Test', dcParams);
    expect(dc.credentials.parameters['login'], dcParams.login);
    expect(dc.credentials.parameters['password'], dcParams.password);
    dc = await grok.dapi.connections.save(dc);
    expect((dc.parameters as any)['schema'], null);
    expect((dc.parameters as any)['db'], dcParams.db);
    expect(dc.friendlyName, 'Local DG Test');
    expect((await grok.dapi.connections.find(dc.id)).id, dc.id);
    // changing credentials
    dc.credentials.parameters['login'] = 'changed_login';
    dc = await grok.dapi.connections.save(dc);
    expect(dc.credentials.openParameters['login'], 'changed_login');

    // changing credentials forEntity
    let credentials = await grok.dapi.credentials.forEntity(dc);
    credentials.parameters['login'] = 'datagrok_dev';
    credentials = await grok.dapi.credentials.save(credentials);
    expect(credentials.openParameters['login'], 'datagrok_dev');

    await grok.dapi.connections.delete(dc);
    expect(await grok.dapi.connections.find(dc.id) == undefined);
  }, {stressTest: true});

  test('JS postprocess', async () => {
    const script = `//language: javascript
//input: dataframe result
//output: dataframe modified
result.columns.remove('id');
modified = result;
`;
    const dc = (await grok.dapi.connections.filter('Datagrok').list())[0];
    const q = dc.query('JS postprocess query test', `SELECT
  i AS id,
  'name_' || i AS name,
  'category_' || (i % 3) AS category,
  i * 10 AS value,
  (i * 1.5)::numeric(10,2) AS price,
  (i % 2 = 0) AS is_active,
  now() - (i || ' days')::interval AS created_at,
  now() AS updated_at,
  'desc_' || i AS description,
  md5(i::text) AS hash
FROM generate_series(1, 10) AS s(i);
    `);
    q.postProcessScript = script;
    q.newId();
    const query = await grok.dapi.queries.save(q);
    expect(query.postProcessScript, script);
    const df = await query.executeTable();
    expect(df.rowCount, 10);
    expect(df.columns.length, 9);
    await grok.dapi.queries.delete(query);
  });

  after(async () => {
    const connections: _DG.DataConnection[] = await grok.dapi.connections.filter(`name="Local DG Test"`).list();
    for (const conn of connections) {
      try {
        await grok.dapi.connections.delete(conn);
      } catch (_) {}
    }
  });
}, { owner: 'ppolovyi@datagrok.ai'});

category('Dapi: connection cache', () => {
  const testFilePath1: string = 'System:AppData/ApiTests/test_files.txt';
  const testFilePath2: string = 'System:AppData/ApiTests/renamed_test_files.txt';

  before(async () => {
    const appData: _DG.DataConnection = await grok.dapi.connections.filter(`shortName="AppData"`).first();
    await grok.functions.call('DropConnectionCache', { 'connection': appData });
  });

  test('Invalidation, performance', async () => {
    // write file to trigger cache bump
    await grok.dapi
      .files.writeAsText(testFilePath1, 'Hello World!');
    // measure first execution time
    let start = Date.now();
    let list = await grok.dapi.files.list('System:AppData/ApiTests');
    const first = Date.now() - start;
    // check if cache was bumped
    expect(list.some((f) => f.name === 'test_files.txt'));
    const second = await getExecutionTime(async () => {
      await grok.dapi.files.list('System:AppData/ApiTests');
    });
    // second execution should be faster
    expect(second * 1.5 < first);

    // cache should be bumped after renaming
    await grok.dapi.files.rename(testFilePath1, 'renamed_test_files.txt');
    list = await grok.dapi.files.list('System:AppData/ApiTests');
    expect(list.some((f) => f.name === 'renamed_test_files.txt'));

    // cache should be bumped after delete
    await grok.dapi.files.delete(testFilePath2);
    list = await grok.dapi.files.list('System:AppData/ApiTests');
    expect(list.every((f) => f.name !== 'renamed_test_files.txt'));
  }, {skipReason: typeof process !== 'undefined' ? 'NodeJS environment' : undefined});

  test('Dataframe: Ids', async () => {
    // not from cache
    const table1 = (await grok.dapi.files.readBinaryDataFrames('System:AppData/ApiTests/datasets/demog.csv'))[0];
    // from cache
    const table2 = (await grok.dapi.files.readBinaryDataFrames('System:AppData/ApiTests/datasets/demog.csv'))[0];
    // id should be absent when we read as csv, and second time from cache
    expect(!table1.id && !table2.id, true);
  }, {skipReason: typeof process !== 'undefined' ? 'NodeJS environment' : undefined });

  test('Sequential stress test', async () => {
    const times = DG.Test.isInBenchmark ? 100 : 10;
    let demogCsvReads1 = [];
    await grok.dapi.files.readCsv('System:AppData/ApiTests/datasets/demog.csv');
    await grok.dapi.files.readCsv('System:AppData/ApiTests/cars.csv');
    for (let i = 0; i < times; i++)
      demogCsvReads1.push(await getExecutionTime(async () => {
        await grok.dapi.files.readCsv('System:AppData/ApiTests/datasets/demog.csv');
      }));

    let carsReads1 = [];
    for (let i = 0; i < times; i++)
      carsReads1.push(await getExecutionTime(async () => {
        await grok.dapi.files.readCsv('System:AppData/ApiTests/cars.csv');
      }));

    let carsReads2 = [];
    let demogCsvReads2 = [];
    for (let i = 0; i < times; i++) {
      demogCsvReads2.push(await getExecutionTime(async () => {
        await grok.dapi.files.readCsv('System:AppData/ApiTests/datasets/demog.csv');
      }));
      carsReads2.push(await getExecutionTime(async () => {
        await grok.dapi.files.readCsv('System:AppData/ApiTests/cars.csv');
      }));
    }
    const demog1Median = median(demogCsvReads1);
    const demog2Median = median(demogCsvReads2);
    expect(demog2Median < demog1Median * 2, true);

    const cars1Median = median(carsReads1);
    const cars2Median = median(carsReads2);
    expect(cars2Median < cars1Median * 5, true);
  }, { benchmark: true, timeout: 120000, skipReason: typeof process !== 'undefined' ? 'NodeJS environment' : undefined });

  after(async () => {
    try {
      await grok.dapi.files.delete(testFilePath1);
    } catch (_) { }
    try {
      await grok.dapi.files.delete(testFilePath2);
    } catch (_) { }
  });
}, { owner: 'ppolovyi@datagrok.ai'});

category('Dapi: TableQuery', () => {
  let dc: _DG.DataConnection;
  const tableName = 'public.orders';
  const fields = ['orderid', 'freight'];
  const whereClauses = [{
    field: 'orderid',
    pattern: '10250',
    dataType: DG.TYPE.STRING
  }];
  const aggregationsDb = [{
    colName: 'orderid',
    aggType: 'count',
  }];
  const havingDb = [{
    field: 'COUNT(shipcountry)',
    pattern: '2',
    dataType: DG.TYPE.STRING
  }];
  const orderByDb = [{
    field: 'orderid',
  }];
  let fromTable: _DG.TableInfo;
  let from: string;

  before(async () => {
    fromTable = DG.TableInfo.fromDataFrame(grok.data.testData('demog', 5000));
    from = fromTable.name;
    const dcParams = {
      dataSource: 'Postgres', server: 'dev.datagrok.ai:54322', db: 'northwind',
      login: 'datagrok', password: 'datagrok'
    };
    dc = DG.DataConnection.create('test', dcParams);
    dc = await grok.dapi.connections.save(dc);
  });

  after(async () => {
    await grok.dapi.connections.delete(dc);
  });

  test('Create', async () => {
    const tq = DG.TableQuery.create(dc);
    expect(tq instanceof DG.TableQuery, true);
  });

  test('Table', async () => {
    const tq = DG.TableQuery.create(dc);
    tq.table = tableName;
    expect(tq.table, tableName);
  });

  test('Fields', async () => {
    const tq = DG.TableQuery.create(dc);
    tq.fields = fields;
    expectArray(tq.fields, fields);
  });

  test('Where clauses', async () => {
    const tq = DG.TableQuery.create(dc);
    tq.fields = fields;
    tq.where = whereClauses;
    expectArray(tq.where, whereClauses);
  });

  test('Aggregations', async () => {
    const tq = DG.TableQuery.create(dc);
    tq.aggregations = aggregationsDb;
    expectArray(tq.aggregations, aggregationsDb);
  });

  test('Having', async () => {
    const tq = DG.TableQuery.create(dc);
    tq.having = havingDb;
    expectArray(tq.having, havingDb);
  });

  test('Order by', async () => {
    const tq = DG.TableQuery.create(dc);
    tq.orderBy = orderByDb;
    expectArray(tq.orderBy, orderByDb);
  });

  test('From table', async () => {
    const dtqb = DG.TableQuery.fromTable(fromTable);
    expect(dtqb instanceof DG.TableQueryBuilder, true);
  });

  test('From', async () => {
    const dtqb = DG.TableQuery.from(from);
    expect(dtqb instanceof DG.TableQueryBuilder, true);
  });
}, { owner: 'ppolovyi@datagrok.ai'});

/*
category('Dapi: TableQueryBuilder', () => {
  before(async () => {
    table = grok.data.testData('demog', 5000);
    fromTable = DG.TableInfo.fromDataFrame(table);
    from = fromTable.name;
  });

  let fromTable: DG.TableInfo;
  let from: string;
  let table: DG.DataFrame;
  const fields = ['race'];

  test('From table', async () => {
    const dtqb = DG.TableQueryBuilder.fromTable(fromTable);
    expect(dtqb instanceof DG.TableQueryBuilder, true);
  });

  test('From', async () => {
    const dtqb = DG.TableQueryBuilder.from(from);
    expect(dtqb instanceof DG.TableQueryBuilder, true);
  });

  test('Select all', async () => {
    let dtqb = DG.TableQueryBuilder.fromTable(fromTable);
    dtqb = dtqb.selectAll();
    const tq = dtqb.build();
    expectArray(tq.fields, table.columns.names());
  });
  
  test('Select', async () => {
    let dtqb = DG.TableQueryBuilder.fromTable(fromTable);
    dtqb = dtqb.select(fields);
    const tq = dtqb.build();
    expectArray(tq.fields, fields);
  });

  test('Group by', async () => {
    let dtqb = DG.TableQueryBuilder.fromTable(fromTable);
    dtqb = dtqb.groupBy(fields);
    dtqb.build();
  });

  test('Pivot on', async () => {
    let dtqb = DG.TableQueryBuilder.fromTable(fromTable);
    dtqb = dtqb.pivotOn(fields);
    dtqb.build();
  });

  test('Where', async () => {
    let dtqb = DG.TableQueryBuilder.fromTable(fromTable);
    dtqb = dtqb.where('race', 'Asian');
    const tq = dtqb.build();
    expectObject(tq.where[0], {field: 'race', pattern: 'Asian'});
  });

  test('Sort by', async () => {
    let dtqb = DG.TableQueryBuilder.fromTable(fromTable);
    dtqb = dtqb.sortBy('age');
    dtqb.build();
  });

  test('Limit', async () => {
    let dtqb = DG.TableQueryBuilder.fromTable(fromTable);
    dtqb = dtqb.limit(10);
    dtqb.build();
  });

  test('Build', async () => {
    const dtqb = DG.TableQueryBuilder.fromTable(fromTable);
    const tq = dtqb.build();
    expect(tq instanceof DG.TableQuery, true);
  });
});
*/

async function getExecutionTime(f: () => any) {
  const start = Date.now();
  await f();
  return Date.now() - start;
}

function median(numbers: number[]) {
  const sorted = Array.from(numbers).sort((a, b) => a - b);
  const middle = Math.floor(sorted.length / 2);
  if (sorted.length % 2 === 0)
    return (sorted[middle - 1] + sorted[middle]) / 2;
  return sorted[middle];
}
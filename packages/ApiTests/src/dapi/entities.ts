import {after, before, category, expect, test, expectArray} from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// category('Dapi: properties', () => {
//
//     test('Dapi: properties - save, get, delete', async () => {
//         let group = DG.Group.create('js-api-test-group1');
//         group = await grok.dapi.groups.save(group);
//
//         let properties = {
//             'entityId': group.id,
//             'property': 'myProp',
//             'value': 'value'
//         };
//
//         await grok.dapi.entities.saveProperties([properties]);
//         await grok.dapi.entities.getProperties(group);
//         await grok.dapi.entities.deleteProperties([properties]);
//     });
//
// });

category('DataConnection', () => {
  
});

category('TableQuery', async () => {
  let dc: DG.DataConnection;
  const tableName = 'public.orders';
  const fields = ['orderid', 'freight'];
  const whereClauses = [{orderid: '10250'}];
  const aggregationsDb = [{orderid: '10250'}];
  const havingDb = [{'COUNT(shipcountry)': '2'}];
  const orderByDb = [{orderid: 'ASC'}];
  const fromTable = await grok.dapi.tables.first();
  const from = fromTable.name;

  before(async () => {
    const dcParams = {dataSource: 'PostgresNet', server: 'dev.datagrok.ai:54322', db: 'northwind',
      login: 'datagrok', password: 'datagrok'};
    dc = DG.DataConnection.createDB('test', dcParams);
    dc = await grok.dapi.connections.save(dc);
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
    tq.whereClauses = whereClauses;
    expectArray(tq.whereClauses, whereClauses);
  });

  test('Aggregations', async () => {
    const tq = DG.TableQuery.create(dc);
    tq.aggregationsDb = aggregationsDb;
    expectArray(tq.aggregationsDb, aggregationsDb);
  });

  test('Having', async () => {
    const tq = DG.TableQuery.create(dc);
    tq.havingDb = havingDb;
    expectArray(tq.havingDb, havingDb);
  });

  test('Order by', async () => {
    const tq = DG.TableQuery.create(dc);
    tq.orderByDb = orderByDb;
    expectArray(tq.orderByDb, orderByDb);
  });

  test('From table', async () => {
    const dtqb = DG.TableQuery.fromTable(fromTable);
    expect(dtqb instanceof DG.DbTableQueryBuilder, true);
  });

  test('From', async () => {
    const dtqb = DG.TableQuery.from(from);
    expect(dtqb instanceof DG.DbTableQueryBuilder, true);
  });
});

category('DbTableQueryBuilder', async () => {
  const fromTable = await grok.dapi.tables.first();
  const from = fromTable.name;

  test('From table', async () => {
    const dtqb = DG.DbTableQueryBuilder.fromTable(fromTable);
    expect(dtqb instanceof DG.DbTableQueryBuilder, true);
  });

  test('From', async () => {
    const dtqb = DG.DbTableQueryBuilder.from(from);
    expect(dtqb instanceof DG.DbTableQueryBuilder, true);
  });

  test('Select all', async () => {
    let dtqb = DG.DbTableQueryBuilder.fromTable(fromTable);
    dtqb.selectAll();
  });
  
  test('Select', async () => {
    let dtqb = DG.DbTableQueryBuilder.fromTable(fromTable);
    dtqb.select([]);
  });

  test('Group by', async () => {
    let dtqb = DG.DbTableQueryBuilder.fromTable(fromTable);
    dtqb.groupBy([]);
  });

  test('Pivot on', async () => {
    let dtqb = DG.DbTableQueryBuilder.fromTable(fromTable);
    dtqb.pivotOn([]);
  });

  test('where', async () => {
    let dtqb = DG.DbTableQueryBuilder.fromTable(fromTable);
    dtqb.where('', '');
  });

  test('Sort by', async () => {
    let dtqb = DG.DbTableQueryBuilder.fromTable(fromTable);
    dtqb.sortBy('');
  });

  test('Limit', async () => {
    let dtqb = DG.DbTableQueryBuilder.fromTable(fromTable);
    dtqb.limit(10);
  });

  test('Build', async () => {

  });
});

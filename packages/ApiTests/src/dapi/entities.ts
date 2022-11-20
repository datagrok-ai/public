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

category('TableQuery', () => {
  let dc: DG.DataConnection;
  const tableName = 'public.orders';
  const fields = ['orderid', 'freight'];
  const whereClauses = [{
    field: 'orderid',
    pattern: '10250',
  }];
  const aggregationsDb = [{
    colName: 'orderid',
    aggType: 'count',
  }];
  const havingDb = [{
    field: 'COUNT(shipcountry)',
    pattern: '2',
  }];
  const orderByDb = [{
    field: 'orderid',
  }];
  let fromTable: DG.TableInfo;
  let from: string;

  before(async () => {
    fromTable = await grok.dapi.tables.first();
    from = fromTable.name;
    const dcParams = {dataSource: 'PostgresNet', server: 'dev.datagrok.ai:54322', db: 'northwind',
      login: 'datagrok', password: 'datagrok'};
    dc = DG.DataConnection.create('test', dcParams);
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
    expect(dtqb instanceof DG.DbTableQueryBuilder, true);
  });

  test('From', async () => {
    const dtqb = DG.TableQuery.from(from);
    expect(dtqb instanceof DG.DbTableQueryBuilder, true);
  });
});

category('DbTableQueryBuilder', () => {
  before(async () => {
    fromTable = await grok.dapi.tables.first();
    from = fromTable.name;
  });

  let fromTable: DG.TableInfo;
  let from: string;

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

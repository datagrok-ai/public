import {after, before, category, expect, test} from "@datagrok-libraries/utils/src/test";
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
  const tableName = 'public.chats';
  const fields = ['id', 'name'];
  const fieldsAll = ['id', 'name', 'direct', 'firendly_name', 'group_id', 'private', 'views'];

  before(async () => {
    const dcParams = {
      dataSource: 'PostgreSQL', server: 'localhost:5432', db: 'datagrok_dev', login: 'datagrok_dev', password: '123'};
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

  });

  test('Aggregations', async () => {

  });

  test('Having', async () => {

  });

  test('Order by', async () => {

  });

  test('From table', async () => {

  });

  test('From', async () => {

  });
});

category('DbTableQueryBuilder', () => {
  test('From table', async () => {

  });

  test('From', async () => {

  });

  test('Select all', async () => {

  });
  
  test('Select', async () => {

  });

  test('Group by', async () => {

  });

  test('Pivot on', async () => {

  });

  test('where', async () => {

  });

  test('Sort by', async () => {

  });

  test('Limit', async () => {

  });

  test('Build', async () => {

  });
});

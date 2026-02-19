import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {DbColumn, DbSchema, DbTable} from "./cruddy/table";
import {DbEntity, DbEntityType} from "./cruddy/entity";
import {DbQueryEntityCrud} from "./cruddy/crud";
import {CruddyConfig} from "./cruddy/app";

const customersTable = new DbTable({
  name: 'customers',
  columns: [
    new DbColumn({name: 'id', type: 'int', isKey: true}),
    new DbColumn({name: 'firstName', type: 'string'}),
    new DbColumn({name: 'lastName', type: 'string'}),
  ]
});

const ordersTable = new DbTable({
  name: 'orders',
  columns: [
    new DbColumn({name: 'id', type: 'int', isKey: true}),
    new DbColumn({name: 'customerId', type: 'int', ref: 'customer.id'}),
    new DbColumn({name: 'details', type: 'string'}),
  ]
});

const northwindConfig = new CruddyConfig({
  connection: 'Samples:PostgresNorthwind',
  schema: new DbSchema('', [ customersTable, ordersTable ])
});


const customerType = new DbEntityType({type: 'customer', table: customersTable});

const john = new DbEntity(customerType, { id: 7, firstName: 'John', lastName: 'Snow'});

category('SQL construction', () => {

  const crud = new DbQueryEntityCrud('conn', customerType);

  const expectSql = (actual: string, expected: string) => expect(actual.replaceAll('\n', ''), expected);

  test('create', async () => {
    expectSql(crud.getCreateSql(john), "insert into customers (id, firstName, lastName) values (7, 'John', 'Snow')");
  });

  test('read', async () => {
    expectSql(crud.getReadSql(), 'select id, firstName, lastName from customers');
  });

  test('update', async () => {
    expectSql(crud.getUpdateSql(john), "update customers set firstName = 'John', lastName = 'Snow' where id = 7");
  });

  test('delete', async () => {
    expectSql(crud.getDeleteSql(john), "delete from customers where id = 7");
  });

})


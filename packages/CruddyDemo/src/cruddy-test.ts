import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import {DbColumn, DbEntity, DbEntityType, DbQueryEntityCrud, DbTable} from "./cruddy";

const customersTable = new DbTable('', 'customers');

const customersColumns = [
  new DbColumn(customersTable, 'id', 'int', true),
  new DbColumn(customersTable, 'firstName', 'string'),
  new DbColumn(customersTable, 'lastName', 'string')
];

const customerType = new DbEntityType('customer', customersColumns);

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
    expectSql(crud.getUpdateSql(john), "update customers set firstName = 'John' and lastName = 'Snow' where id = 7");
  });

  test('delete', async () => {
    expectSql(crud.getDeleteSql(john), "delete from customers where id = 7");
  });

})


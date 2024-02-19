import * as cruddy from './cruddy';
import {DbColumn, DbEntityType, DbQueryEntityCrud, DbTable} from "./cruddy";

const customersTable = new DbTable('', 'customers');

const customersColumns = [
  new DbColumn(customersTable, 'id', 'int'),
  new DbColumn(customersTable, 'firstName', 'string'),
  new DbColumn(customersTable, 'lastName', 'string')
];

const customerType = new DbEntityType('customer', customersColumns);

test('read sql', () => {
  const crud = new DbQueryEntityCrud('conn', customerType);
  expect(crud.getReadSql()).toEqual('select id, firstName, lastName from customers');
});
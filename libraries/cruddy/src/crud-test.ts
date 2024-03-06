import * as cruddy from './cruddy';
import {DbColumn, DbEntityType, DbQueryEntityCrud, DbTable} from "./cruddy";

const customersTable = new DbTable('', 'customers');

const customersColumns = [
  new DbColumn(customersTable, 'id', 'int'),
  new DbColumn(customersTable, 'firstName', 'string'),
  new DbColumn(customersTable, 'lastName', 'string')
];

const customerType = new DbEntityType('customer', customersColumns);

const crud = new DbQueryEntityCrud('conn', customerType);
console.log(crud.getReadSql());

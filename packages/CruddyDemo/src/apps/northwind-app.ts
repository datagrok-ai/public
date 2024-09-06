import {DbColumn, DbSchema, DbTable} from "../cruddy/table";
import {DbEntityType} from "../cruddy/entity";
import {CruddyApp, CruddyConfig, CruddyEntityView} from "../cruddy/app";


const customersTable = new DbTable({
  name: 'customers',
  columns: [
    new DbColumn({name: 'customerid', type: 'string', isKey: true}),
    new DbColumn({name: 'companyname', type: 'string'}),
    new DbColumn({name: 'city', type: 'string'}),
    new DbColumn({name: 'region', type: 'string', nullable: true}),
    new DbColumn({name: 'fax', type: 'string'}),
    new DbColumn({name: 'contacttitle', type: 'string'}),
    new DbColumn({name: 'country', type: 'string'}),
    new DbColumn({name: 'postalcode', type: 'string'}),
    new DbColumn({name: 'contactname', type: 'string'}),
    new DbColumn({name: 'address', type: 'string'}),
    new DbColumn({name: 'phone', type: 'string', nullable: true}),
  ]
});

const employeesTable = new DbTable({
  name: 'employees',
  columns: [
    new DbColumn({name: 'employeeid', type: 'int', isKey: true}),
    new DbColumn({name: 'city', type: 'string'}),
    new DbColumn({name: 'country', type: 'string'}),
    new DbColumn({name: 'title', type: 'string'}),
    new DbColumn({name: 'firstname', type: 'string'}),
    new DbColumn({name: 'lastname', type: 'string'}),
  ]
});

const categoriesTable = new DbTable({
  name: 'categories',
  columns: [
    new DbColumn({name: 'categoryid', type: 'int', isKey: true}),
    new DbColumn({name: 'categoryname', type: 'string'}),
  ]
});

const productsTable = new DbTable({
  name: 'products',
  columns: [
    new DbColumn({name: 'productid', type: 'int', isKey: true}),
    new DbColumn({name: 'productname', type: 'string'}),
    new DbColumn({name: 'categoryid', type: 'int', ref: 'categories.categoryid'}),
    new DbColumn({name: 'unitprice', type: 'double'}),
  ]
});

const ordersTable = new DbTable({
  name: 'orders',
  columns: [
    new DbColumn({name: 'orderid', type: 'int', isKey: true}),
    new DbColumn({name: 'shippeddate', type: 'datetime'}),
    new DbColumn({name: 'shipcountry', type: 'string'}),
    new DbColumn({name: 'shipregion', type: 'string', nullable: true}),
    new DbColumn({name: 'shipcity', type: 'string'}),
    new DbColumn({name: 'employeeid', type: 'int', ref: 'employees.employeeid'}),
    new DbColumn({name: 'customerid', type: 'string', ref: 'customers.customerid'}),
  ]
});

const orderDetailsTable = new DbTable({
  name: 'order_details',
  columns: [
    new DbColumn({name: 'orderid', type: 'int', ref: 'orders.orderid'}),
    new DbColumn({name: 'productid', type: 'int', ref: 'products.productid'}),
    new DbColumn({name: 'unitprice', type: 'double'}),
    new DbColumn({name: 'quantity', type: 'int'}),
  ]
});

const northwindSchema: DbSchema = new DbSchema('northwind',
  [customersTable, employeesTable, categoriesTable, productsTable, ordersTable, orderDetailsTable]);

export const northwindConfig = new CruddyConfig({
  connection: 'Samples:PostgresNorthwind',
  schema: northwindSchema,
  entityTypes: [
    new DbEntityType({ type: 'Category', table: categoriesTable }),
    new DbEntityType({ type: 'Product', table: productsTable }),
    new DbEntityType({ type: 'Customer', table: customersTable }),
    new DbEntityType({ type: 'Employee', table: employeesTable }),
    new DbEntityType({ type: 'Order', table: ordersTable,
      gridColumnsNames: [
        'orderid', 'shippeddate', 'shipcountry', 'shipregion', 'shipcity',
        'employees.firstname', 'customers.companyname'
      ],
      filters: [
        { type: 'distinct', column: 'shipcountry'},
        { type: 'combo', column: 'shipcity'},
        { type: 'range', column: 'orderid'},
        { type: 'expression', column: 'shipcity'},
    ]}),
    new DbEntityType({ type: 'Order Details', table: orderDetailsTable }),
  ]
});


export const northwindApp = new CruddyApp(
  northwindConfig, [
    //new CruddyEntityView()
  ]
);
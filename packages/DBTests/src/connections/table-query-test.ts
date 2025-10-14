import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';

category('TableQueryBuilder', () => {
    const testTableName: string = 'products';
    let northwind: DG.DataConnection;

    before(async () => {
        northwind = await grok.functions.eval('DbTests:PostgresTest');
    });

    test('no fields', async () => {
        const tq = northwind.buildQuery(testTableName).build();
        const df: DG.DataFrame = await tq.executeTable();
        expect(df.rowCount, 77);
        expect(df.columns.length, 10);
    });

    test('fields limit', async () => {
        const tq = northwind.buildQuery(testTableName)
            .select(['productname', 'unitprice'])
            .limit(10)
            .build();
        const df: DG.DataFrame = await tq.executeTable();
        expect(df.rowCount, 10);
        expect(df.columns.length, 2);
    });

    test('fields where double pattern limit', async () => {
        const tq = northwind.buildQuery(testTableName)
            .select(['productname', 'unitprice'])
            .where('unitprice', '<=13',  DG.TYPE.FLOAT)
            .limit(3)
            .build();
        const df: DG.DataFrame = await tq.executeTable();
        expect(df.rowCount, 3);
        expect(df.columns.length, 2);
    });

    test('fields where datetime pattern', async () => {
        const tq = northwind.buildQuery('orders')
            .select(['shipcity', 'shipcountry', 'orderdate'])
            .where('orderdate', 'before 1996-08-26',  DG.TYPE.DATE_TIME)
            .build();
        const df: DG.DataFrame = await tq.executeTable();
        expect(df.rowCount, 42);
        expect(df.columns.length, 3);
    });

    test('fields where string pattern', async () => {
        const tq = northwind.buildQuery('customers')
            .where('region', 'is not null')
            .build();
        const df: DG.DataFrame = await tq.executeTable();
        expect(df.rowCount, 31);
        expect(df.columns.length, 11);
    });

    test('aggregation max sort', async () => {
        // translates to select * from products
        const tq = northwind.buildQuery(testTableName)
            .max('unitprice', 'max price')
            .groupBy(['categoryid', 'unitprice'])
            .sortBy('unitprice', false)
            .build();
        const df: DG.DataFrame = await tq.executeTable();
        expect(df.rowCount, 73);
        expect(df.columns.length, 3);
        expect(df.getCol('max price').get(0), 263.50);
    });

    test('aggregation max sort left join', async () => {
        // translates to select * from products
        const tq = northwind.buildQuery(testTableName)
            .leftJoin('categories', ['categoryid'], ['categoryid'])
            .max('unitprice', 'max price')
            .groupBy(['categoryname', 'unitprice'])
            .sortBy('unitprice', false)
            .build();
        const df: DG.DataFrame = await tq.executeTable();
        expect(df.rowCount, 73);
        expect(df.columns.length, 3);
        expect(df.getCol('categoryname').get(0), 'Beverages');
    });

    test('aggregation max sort left join having', async () => {
        // translates to select * from products
        const tq = northwind.buildQuery(testTableName)
            .leftJoin('categories', ['categoryid'], ['categoryid'])
            .max('unitprice', 'max price')
            .groupBy(['categoryname', 'unitprice'])
            .having('unitprice', '>100', DG.TYPE.FLOAT)
            .sortBy('unitprice', false)
            .build();
        const df: DG.DataFrame = await tq.executeTable();
        expect(df.rowCount, 2);
        expect(df.columns.length, 3);
        expect(df.getCol('categoryname').get(0), 'Beverages');
    });

    test('aggregation max left inner joins having', async () => {
        // translates to select * from products
        const tq = northwind.buildQuery(testTableName)
            .leftJoin('categories', ['categoryid'], ['categoryid'])
            .innerJoin('suppliers', ['supplierid'], ['supplierid'])
            .valueCount('categoryname', 'categories by supplier')
            .groupBy(['companyname', 'products.supplierid', 'country'])
            .having('country', 'in (USA, UK)')
            .build();
        const df: DG.DataFrame = await tq.executeTable();
        expect(df.rowCount, 6);
        expect(df.columns.length, 4);
    });
});

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

category('Database Meta: DbInfo', async () => {
    let conn: DG.DataConnection;
    let dbInfo: DG.DbInfo;
    const testComment = 'test-comment-' + Date.now();
    const testLlmComment = 'llm-comment-' + Date.now();

    before(async () => {
        conn = await grok.functions.eval('DbTests:PostgreSQLDBTests');
        dbInfo = (await grok.data.db.getInfo(conn, 'public'))[0];
        try {
            await dbInfo?.clearProperties();
        } catch (_) {}
    });

    test('Base properties', async () => {
        expect(dbInfo.connection.id, conn.id);
        expect(dbInfo.comment === undefined || dbInfo.comment === null, true);
        await dbInfo.setComment(testComment);
        expect(dbInfo.comment, testComment);
        let updated = (await grok.data.db.getInfo(conn, 'public'))[0];
        expect(updated.comment, testComment);

        await dbInfo.setLlmComment(testLlmComment);
        expect(dbInfo.llmComment, testLlmComment);
        updated = (await grok.data.db.getInfo(conn, 'public'))[0];
        expect(updated.llmComment, testLlmComment);
    });

    test('Add relation 1', async () => {
        const fromTable: string = 'orders';
        const fromColumns: string[] = ['customer_id'];
        const toTable: string = 'customers';
        const toColumns: string[] = ['id'];
        const props: DG.DbRelationProperties = {
            comment: 'test',
            llmComment: 'test llm comment',
            cardinality: 'one-to-many',
            fromSchema: 'sales',
            toSchema: 'sales'
        };

        const rel = await dbInfo.addRelation(fromTable, fromColumns, toTable, toColumns, props);
        expect(rel.comment, props.comment);
        expect(rel.llmComment, props.llmComment);
        expect(rel.cardinality, props.cardinality);
        expect(rel.fromTable, fromTable);
        expect(rel.fromSchema, props.fromSchema);
        expect(rel.fromColumns, fromColumns);
        expect(rel.toTable, toTable);
        expect(rel.toSchema, props.toSchema);
        expect(rel.toColumns, toColumns);

        let all = await dbInfo.getRelations();
        const found = all.find(r => r.fromTable === fromTable);
        expect(found !== null, true);
    });

    test('Add relation 2', async () => {
        const fromTable: string = 'orders';
        const fromColumns: string[] = ['customer_name', 'customer_address'];
        const toTable: string = 'customers';
        const toColumns: string[] = ['name', 'address'];

        const props: DG.DbRelationProperties = {
            comment: 'test',
            llmComment: 'test llm comment',
            cardinality: 'one-to-many',
            fromSchema: 'sales',
            toSchema: 'sales'
        };

        const rel = await dbInfo.addRelation(fromTable, fromColumns, toTable, toColumns, props);
        expect(rel.fromColumns, fromColumns);
        expect(rel.toColumns, toColumns);

        let all = await dbInfo.getRelations();
        const found = all.find(r => r.fromTable === fromTable);
        expect(found !== null, true);
    });

    test('Clear properties', async () => {
        await dbInfo.setComment(testComment);
        const fromTable: string = 'orders';
        const fromColumns: string[] = ['customer_name'];
        const toTable: string = 'customers';
        const toColumns: string[] = ['name'];
        const props: DG.DbRelationProperties = {
            comment: 'test',
            llmComment: 'test llm comment',
            cardinality: 'one-to-many'
        };

        await dbInfo.addRelation(fromTable, fromColumns, toTable, toColumns, props);
        await dbInfo.clearProperties();
        let updated = (await grok.data.db.getInfo(conn, 'public'))[0];
        expect(updated.comment == null, true);
        const rels = await dbInfo.getRelations();
        expect(rels.length, 0);
    });

    after(async () => {
        try {
            await dbInfo?.clearProperties();
        } catch (_) {}
    });
});

category('Database Meta: DbSchemaInfo', async () => {
    let conn: DG.DataConnection;
    let dbInfo: DG.DbInfo;
    const testComment = 'test-schema-comment-' + Date.now();
    const testLlmComment = 'llm-schema-comment-' + Date.now();

    before(async () => {
        conn = await grok.functions.eval('DbTests:PostgreSQLDBTests');
        dbInfo = (await grok.data.db.getInfo(conn, 'public'))[0];
        try {
            await dbInfo?.clearProperties();
        } catch (_) {}
    });

    test('Base properties', async () => {
        const schemas = await dbInfo.getSchemas();
        expect(schemas.length, 3);
        let publicSchema = schemas.find((s) => s.name === 'public');
        expect(publicSchema !== null, true);
        publicSchema = publicSchema!;
        expect(publicSchema.connection.id, conn.id);

        await publicSchema.setComment(testComment);
        expect(publicSchema.comment, testComment);
        let updated = (await grok.data.db.getInfo(conn, 'public'))[0];
        let updatedSchema = (await updated.getSchemas()).find((s) => s.name === 'public')!;
        expect(updatedSchema.comment, testComment);

        await publicSchema.setLlmComment(testLlmComment);
        expect(publicSchema.llmComment, testLlmComment);
        updated = (await grok.data.db.getInfo(conn, 'public'))[0];
        updatedSchema = (await updated.getSchemas()).find((s) => s.name === 'public')!;
        expect(updatedSchema.llmComment, testLlmComment);
    });

    test('Annotate table', async () => {
        const schemas = await dbInfo.getSchemas();
        let publicSchema = schemas.find((s) => s.name === 'public')!;
        const tables = await publicSchema.getTables();
        const t = tables.find((t) => t.name == 'mock_data') ?? tables[0];

        const props: DG.DbTableProperties = {
            comment: 'table-comment',
            rowCount: 123,
        };

        await publicSchema.annotateTable(t, props);

        const refreshed = await dbInfo.getSchemas();
        const s2 = refreshed.find(x => x.name === publicSchema.name)!;
        const t2 = (await s2.getTables()).find(x => x.name === t.name)!;
        expect(t2.tags[DG.Tags.DbComment], 'table-comment');
        expect(t2.tags[DG.Tags.DbTableRowCount], 123);

        // full refresh
        const updatedInfo = (await grok.data.db.getInfo(conn, 'public'))[0];
        const updatedSchemas = await updatedInfo.getSchemas();
        const s3 = updatedSchemas.find(x => x.name === publicSchema.name)!;
        const t3 = (await s3.getTables()).find(x => x.name === t.name)!;
        expect(t3.tags[DG.Tags.DbComment], 'table-comment');
        expect(t3.tags[DG.Tags.DbTableRowCount], 123);
    });

    test('Annotate column', async () => {
        const schemas = await dbInfo.getSchemas();
        const s = schemas[0];
        const tables = await s.getTables();
        const t = tables[0];
        const col = t.columns[0];

        const props: DG.DbColumnProperties = {
            comment: 'col-comment',
            min: 1,
            max: 10,
        };

        await s.annotateColumn(t, col, props);

        const refreshed = await dbInfo.getSchemas();
        const s2 = refreshed.find(x => x.name === s.name)!;
        const t2 = (await s2.getTables()).find(x => x.name === t.name)!;
        const c2 = t2.columns.find(x => x.name === col.name)!;
        expect(c2.tags[DG.Tags.DbColumnMin], props.min);
        expect(c2.tags[DG.Tags.DbColumnMax], props.max);
        expect(c2.tags[DG.Tags.DbComment], props.comment);

        // full refresh
        const updatedInfo = (await grok.data.db.getInfo(conn, 'public'))[0];
        const updatedSchemas = await updatedInfo.getSchemas();
        const s3 = updatedSchemas.find(x => x.name === x.name)!;
        const t3 = (await s3.getTables()).find(x => x.name === t.name)!;
        const c3 = t3.columns.find(x => x.name === col.name)!;
        expect(c3.tags[DG.Tags.DbColumnMin], props.min);
        expect(c3.tags[DG.Tags.DbColumnMax], props.max);
        expect(c3.tags[DG.Tags.DbComment], props.comment);
    });

    after(async () => {
        try {
            await dbInfo?.clearProperties();
        } catch (_) {}
    });
});

category('Database Meta: DbRelationInfo', async () => {
    let conn: DG.DataConnection;
    let dbInfo: DG.DbInfo;

    before(async () => {
        conn = await grok.functions.eval('DbTests:PostgreSQLDBTests');
        dbInfo = (await grok.data.db.getInfo(conn, 'public'))[0];
        try {
            await dbInfo?.clearProperties();
        } catch (_) {}

        const fromTable: string = 'orders';
        const fromColumns: string[] = ['customer_id'];
        const toTable: string = 'customers';
        const toColumns: string[] = ['id'];

        const props: DG.DbRelationProperties = {
            comment: 'test',
            llmComment: 'test llm comment',
            cardinality: 'one-to-many',
            fromSchema: 'sales',
            toSchema: 'sales',
            isPrimaryPath: true
        };

        await dbInfo.addRelation(fromTable, fromColumns, toTable, toColumns, props);
    });

    test('Base properties', async () => {
        const relations = await dbInfo.getRelations();
        expect(relations.length, 1);
        const rel = relations[0];
        expect(rel.fromTable, 'orders');
        expect(rel.comment, 'test');
        expect(rel.llmComment, 'test llm comment');
        expect(rel.cardinality, 'one-to-many');
        expect(rel.isPrimaryPath, true);
    });

    after(async () => {
        try {
            await dbInfo?.clearProperties();
        } catch (_) {}
    });
});

import {after, before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

category('Database Meta: DbInfo', async () => {
    let conn: DG.DataConnection;
    let dbInfo: DG.DbInfo;
    const testComment = 'test-comment-' + Date.now();
    const testLlmComment = 'llm-comment-' + Date.now();
    const relName1 = 'rel_test1_' + Date.now();
    const relName2 = 'rel_test2_' + Date.now();

    before(async () => {
        conn = await grok.functions.eval('DbTests:PostgreSQLDBTests');
        dbInfo = await grok.dapi.connections.getDatabaseInfo(conn);
        try {
            await dbInfo?.clearProperties();
        } catch (_) {}
    });

    test('Base properties', async () => {
        expect(dbInfo.connection.id, conn.id);
        expect(dbInfo.comment === null, true);
        await dbInfo.setComment(testComment);
        expect(dbInfo.comment, testComment);
        let updated = await grok.dapi.connections.getDatabaseInfo(conn);
        expect(updated.comment, testComment);

        await dbInfo.setLlmComment(testLlmComment);
        expect(dbInfo.llmComment, testLlmComment);
        updated = await grok.dapi.connections.getDatabaseInfo(conn);
        expect(updated.llmComment, testLlmComment);
    });

    test('Add relation 1', async () => {
        const props: DG.DbRelationProperties = {
            comment: 'test',
            llmComment: 'test llm comment',
            cardinality: 'one-to-many',
            fromTable: 'orders',
            fromSchema: 'sales',
            fromColumns: ['customer_id'],
            toTable: 'customers',
            toSchema: 'sales',
            toColumns: ['id'],
        };

        const rel = await dbInfo.addRelation(relName1, props);
        expect(rel.comment, props.comment);
        expect(rel.llmComment, props.llmComment);
        expect(rel.cardinality, props.cardinality);
        expect(rel.fromTable, props.fromTable);
        expect(rel.fromSchema, props.fromSchema);
        expect(rel.fromColumns, props.fromColumns);
        expect(rel.toTable, props.toTable);
        expect(rel.toSchema, props.toSchema);
        expect(rel.toColumns, props.toColumns);

        let all = await dbInfo.relations;
        const found = all.find(r => r.name === relName1);
        expect(found !== null, true);
    });

    test('Add relation 2', async () => {
        const props: DG.DbRelationProperties = {
            comment: 'test',
            llmComment: 'test llm comment',
            cardinality: 'one-to-many',
            fromTable: 'orders',
            fromSchema: 'sales',
            fromColumns: ['customer_name', 'customer_address'],
            toTable: 'customers',
            toSchema: 'sales',
            toColumns: ['name', 'address'],
        };

        const rel = await dbInfo.addRelation(relName2, props);
        expect(rel.fromColumns, props.fromColumns);
        expect(rel.toColumns, props.toColumns);

        let all = await dbInfo.relations;
        const found = all.find(r => r.name === relName2);
        expect(found !== null, true);
    });

    test('Clear properties', async () => {
        await dbInfo.setComment(testComment);
        const props: DG.DbRelationProperties = {
            comment: 'test',
            llmComment: 'test llm comment',
            cardinality: 'one-to-many',
            fromTable: 'orders',
            fromColumns: ['customer_name'],
            toTable: 'customers',
            toColumns: ['name'],
        };

        await dbInfo.addRelation(relName1, props);
        await dbInfo.clearProperties();
        let updated = await grok.dapi.connections.getDatabaseInfo(conn);
        expect(updated.comment == null, true);
        const rels = await dbInfo.relations;
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
        dbInfo = await grok.dapi.connections.getDatabaseInfo(conn);
        try {
            await dbInfo?.clearProperties();
        } catch (_) {}
    });

    test('Base properties', async () => {
        const schemas = await dbInfo.schemas;
        expect(schemas.length, 3);
        let publicSchema = schemas.find((s) => s.name === 'public');
        expect(publicSchema !== null, true);
        publicSchema = publicSchema!;
        expect(publicSchema.connection.id, conn.id);

        await publicSchema.setComment(testComment);
        expect(publicSchema.comment, testComment);
        let updated = await grok.dapi.connections.getDatabaseInfo(conn);
        let updatedSchema = (await updated.schemas).find((s) => s.name === 'public')!;
        expect(updatedSchema.comment, testComment);

        await publicSchema.setLlmComment(testLlmComment);
        expect(publicSchema.llmComment, testLlmComment);
        updated = await grok.dapi.connections.getDatabaseInfo(conn);
        updatedSchema = (await updated.schemas).find((s) => s.name === 'public')!;
        expect(updatedSchema.llmComment, testLlmComment);
    });

    test('Annotate table', async () => {
        const schemas = await dbInfo.schemas;
        let publicSchema = schemas.find((s) => s.name === 'public')!;
        const tables = await publicSchema.tables;
        const t = tables.find((t) => t.name == 'mock_data') ?? tables[0];

        const props: DG.DbTableProperties = {
            comment: 'table-comment',
            rowCount: 123,
        };

        await publicSchema.annotateTable(t, props);

        const refreshed = await dbInfo.schemas;
        const s2 = refreshed.find(x => x.name === publicSchema.name)!;
        const t2 = (await s2.tables).find(x => x.name === t.name)!;
        expect(t2.tags['comment'], 'table-comment');
        expect(t2.tags['rowCount'], 123);

        // full refresh
        const updatedInfo = await grok.dapi.connections.getDatabaseInfo(conn);
        const updatedSchemas = await updatedInfo.schemas;
        const s3 = updatedSchemas.find(x => x.name === publicSchema.name)!;
        const t3 = (await s3.tables).find(x => x.name === t.name)!;
        expect(t3.tags['comment'], 'table-comment');
        expect(t3.tags['rowCount'], 123);
    });

    test('Annotate column', async () => {
        const schemas = await dbInfo.schemas;
        const s = schemas[0];
        const tables = await s.tables;
        const t = tables[0];
        const col = t.columns[0];

        const props: DG.DbColumnProperties = {
            comment: 'col-comment',
            min: 1,
            max: 10,
        };

        await s.annotateColumn(t, col, props);

        const refreshed = await dbInfo.schemas;
        const s2 = refreshed.find(x => x.name === s.name)!;
        const t2 = (await s2.tables).find(x => x.name === t.name)!;
        const c2 = t2.columns.find(x => x.name === col.name)!;
        expect(c2.tags['min'], props.min);
        expect(c2.tags['max'], props.max);
        expect(c2.tags['comment'], props.comment);

        // full refresh
        const updatedInfo = await grok.dapi.connections.getDatabaseInfo(conn);
        const updatedSchemas = await updatedInfo.schemas;
        const s3 = updatedSchemas.find(x => x.name === x.name)!;
        const t3 = (await s3.tables).find(x => x.name === t.name)!;
        const c3 = t3.columns.find(x => x.name === col.name)!;
        expect(c3.tags['min'], props.min);
        expect(c3.tags['max'], props.max);
        expect(c3.tags['comment'], props.comment);
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
    const relName = 'test_rel_fk';

    before(async () => {
        conn = await grok.functions.eval('DbTests:PostgreSQLDBTests');
        dbInfo = await grok.dapi.connections.getDatabaseInfo(conn);
        try {
            await dbInfo?.clearProperties();
        } catch (_) {}

        const props: DG.DbRelationProperties = {
            comment: 'test',
            llmComment: 'test llm comment',
            cardinality: 'one-to-many',
            fromTable: 'orders',
            fromSchema: 'sales',
            fromColumns: ['customer_id'],
            toTable: 'customers',
            toSchema: 'sales',
            toColumns: ['id'],
            isPrimaryPath: true
        };

        await dbInfo.addRelation(relName, props);
    });

    test('Base properties', async () => {
        const relations = await dbInfo.relations;
        expect(relations.length, 1);
        const rel = relations[0];
        expect(rel.name, relName);

        const unique = 'rel-comment-' + Date.now();
        await rel.setComment(unique);
        await rel.setLlmComment('llm-' + unique);
        await rel.setCardinality('one-to-one');
        await rel.setFromTable('new_orders');
        await rel.setFromSchema('new_sales');
        await rel.setFromColumns(['new_customer_id']);
        await rel.setToTable('new_customer');
        await rel.setFromSchema('new_sales');
        await rel.setFromColumns(['new_id']);
        await rel.setIsPrimaryPath(false);

        let updated = (await dbInfo.relations).find(x => x.name === rel.name)!;
        expect(updated.comment, unique);
        expect(updated.llmComment, 'llm-' + unique);
        expect(updated.cardinality, 'one-to-one');
        expect(updated.fromTable, 'new_orders');
        expect(updated.fromSchema, 'new_sales');
        expect(updated.fromColumns[0], 'new_id');
        expect(updated.isPrimaryPath, false);

        // full refresh
        const updatedInfo = await grok.dapi.connections.getDatabaseInfo(conn);
        updated = (await updatedInfo.relations).find(x => x.name === rel.name)!;
        expect(updated.comment, unique);
        expect(updated.llmComment, 'llm-' + unique);
        expect(updated.cardinality, 'one-to-one');
        expect(updated.fromTable, 'new_orders');
        expect(updated.fromSchema, 'new_sales');
        expect(updated.fromColumns[0], 'new_id');
        expect(updated.isPrimaryPath, false);
    });

    after(async () => {
        try {
            await dbInfo?.clearProperties();
        } catch (_) {}
    });
});

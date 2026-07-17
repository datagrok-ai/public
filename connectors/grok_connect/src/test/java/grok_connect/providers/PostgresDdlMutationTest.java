package grok_connect.providers;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.Provider;
import grok_connect.table_mutation.AddKey;
import grok_connect.table_mutation.AlterTable;
import grok_connect.table_mutation.ColumnSpec;
import grok_connect.table_mutation.CreateIndex;
import grok_connect.table_mutation.CreateTable;
import grok_connect.table_mutation.DeleteRows;
import grok_connect.table_mutation.DropIndex;
import grok_connect.table_mutation.DropTable;
import grok_connect.table_mutation.IndexSpec;
import grok_connect.table_mutation.InsertRows;
import grok_connect.table_mutation.MutationBatch;
import grok_connect.table_mutation.MutationConfirmationRequiredException;
import grok_connect.table_mutation.MutationResult;
import grok_connect.table_mutation.MutationRunner;
import grok_connect.table_mutation.MutationValidationException;
import grok_connect.table_mutation.TableMutation;
import grok_connect.table_mutation.TruncateTable;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

/**
 * DDL execution tests (connector-writes WO-B7) against a containerized Postgres: the seven-op
 * lifecycle with information_schema verification, the dryRun / confirmDestructive / recomputed-plan
 * contract (ARCHITECTURE §3.4.2), live-data pre-checks and hard refusals, transactional-DDL
 * rollback honesty, the single-op-per-request rule, and the capability gate.
 */
class PostgresDdlMutationTest extends ContainerizedProviderBaseTest {
    protected PostgresDdlMutationTest() {
        super(Provider.POSTGRESQL);
    }

    /** Direct JDBC on a SECOND connection, bypassing grok_connect — verification must not depend on the code under test. */
    private void execDirect(String sql) throws SQLException {
        try (Connection c = DriverManager.getConnection(container.getJdbcUrl(), container.getUsername(), container.getPassword());
             Statement statement = c.createStatement()) {
            statement.execute(sql);
        }
    }

    private List<Object[]> queryDirect(String sql) throws SQLException {
        try (Connection c = DriverManager.getConnection(container.getJdbcUrl(), container.getUsername(), container.getPassword());
             PreparedStatement statement = c.prepareStatement(sql);
             ResultSet rs = statement.executeQuery()) {
            List<Object[]> rows = new ArrayList<>();
            int columnCount = rs.getMetaData().getColumnCount();
            while (rs.next()) {
                Object[] row = new Object[columnCount];
                for (int i = 0; i < columnCount; i++)
                    row[i] = rs.getObject(i + 1);
                rows.add(row);
            }
            return rows;
        }
    }

    private long countDirect(String table) throws SQLException {
        return ((Number) queryDirect("SELECT count(*) FROM " + table).get(0)[0]).longValue();
    }

    /** {data_type, is_nullable} from information_schema, or null when the column does not exist. */
    private Object[] columnInfo(String table, String column) throws SQLException {
        List<Object[]> rows = queryDirect("SELECT data_type, is_nullable FROM information_schema.columns"
                + " WHERE table_name = '" + table + "' AND column_name = '" + column + "'");
        return rows.isEmpty() ? null : rows.get(0);
    }

    private boolean tableExists(String table) throws SQLException {
        return !queryDirect("SELECT 1 FROM information_schema.tables WHERE table_name = '" + table + "'").isEmpty();
    }

    private boolean indexExists(String table, String index) throws SQLException {
        return !queryDirect("SELECT 1 FROM pg_indexes WHERE tablename = '" + table
                + "' AND indexname = '" + index + "'").isEmpty();
    }

    @AfterAll
    public void dropScratchTables() throws SQLException {
        // ddl_child first: it may hold an FK into ddl_life
        for (String table : new String[] {"ddl_child", "ddl_life", "ddl_dry_drop", "ddl_confirm", "ddl_nn",
                "ddl_req", "ddl_req_empty", "ddl_ct", "ddl_tx", "ddl_dml"})
            execDirect("DROP TABLE IF EXISTS " + table);
    }

    private MutationResult runMutation(TableMutation mutation) throws Exception {
        mutation.connection = connection;
        FuncCall call = new FuncCall();
        call.func = mutation;
        call.options = new HashMap<>();
        return MutationRunner.execute(provider, call);
    }

    private ColumnSpec spec(String name, String type, boolean nullable, String defaultValue) {
        ColumnSpec c = new ColumnSpec();
        c.name = name;
        c.type = type;
        c.nullable = nullable;
        c.defaultValue = defaultValue;
        return c;
    }

    private AlterTable alter(String table, String action, String columnName) {
        AlterTable m = new AlterTable();
        m.tableName = table;
        m.action = action;
        m.columnName = columnName;
        return m;
    }

    private DropTable dropTable(String table, boolean confirm) {
        DropTable m = new DropTable();
        m.tableName = table;
        m.confirmDestructive = confirm;
        return m;
    }

    private InsertRows insertRows(String table, List<String> columns, List<String> columnTypes, List<List<Object>> rows) {
        InsertRows m = new InsertRows();
        m.tableName = table;
        m.columns = columns;
        m.columnTypes = columnTypes;
        m.rows = rows;
        return m;
    }

    @DisplayName("Lifecycle: createTable(PK+unique index) -> alter ops -> createIndex/addKey(FK)/dropIndex -> truncate -> dropTable, information_schema-verified")
    @Test
    public void lifecycle_createToDrop() throws Exception {
        CreateTable create = new CreateTable();
        create.tableName = "ddl_life";
        create.columns = Arrays.asList(
                spec("id", "int", false, null),
                spec("name", "string", true, null),
                spec("qty", "int", true, "0"));
        create.primaryKey = Collections.singletonList("id");
        create.indexes = Collections.singletonList(new IndexSpec());
        create.indexes.get(0).name = "ux_ddl_life_name";
        create.indexes.get(0).unique = true;
        create.indexes.get(0).columns = Collections.singletonList("name");
        MutationResult result = runMutation(create);
        Assertions.assertNull(result.errorMessage);
        Assertions.assertEquals(0, result.affectedRows);
        Assertions.assertNotNull(result.plan); // audit value on executed DDL too
        Assertions.assertEquals(2, result.plan.statements.size());
        Assertions.assertEquals(Boolean.TRUE, result.plan.transactionalDdl);
        Assertions.assertArrayEquals(new Object[] {"integer", "NO"}, columnInfo("ddl_life", "id"));
        Assertions.assertArrayEquals(new Object[] {"text", "YES"}, columnInfo("ddl_life", "name"));
        Assertions.assertArrayEquals(new Object[] {"integer", "YES"}, columnInfo("ddl_life", "qty"));
        Assertions.assertEquals("0", queryDirect("SELECT column_default FROM information_schema.columns"
                + " WHERE table_name = 'ddl_life' AND column_name = 'qty'").get(0)[0]);
        Assertions.assertTrue(indexExists("ddl_life", "ux_ddl_life_name"));
        execDirect("INSERT INTO ddl_life VALUES (1, 'a', 10), (2, 'b', 20)");

        AlterTable addColumn = alter("ddl_life", "addColumn", null);
        addColumn.column = spec("note", "string", true, null);
        Assertions.assertNull(runMutation(addColumn).errorMessage);
        Assertions.assertArrayEquals(new Object[] {"text", "YES"}, columnInfo("ddl_life", "note"));

        AlterTable rename = alter("ddl_life", "renameColumn", "qty");
        rename.newName = "quantity";
        Assertions.assertNull(runMutation(rename).errorMessage);
        Assertions.assertNull(columnInfo("ddl_life", "qty"));
        Assertions.assertNotNull(columnInfo("ddl_life", "quantity"));

        AlterTable changeType = alter("ddl_life", "changeType", "quantity");
        changeType.newType = "bigint";
        changeType.confirmDestructive = true; // type-change is destructive: pre-checked, confirm required
        result = runMutation(changeType);
        Assertions.assertNull(result.errorMessage);
        Assertions.assertEquals("type-change", result.plan.destructive.get(0).code);
        Assertions.assertEquals("columns.quantity", result.plan.destructive.get(0).path);
        Assertions.assertEquals(2L, result.plan.destructive.get(0).count); // live non-null values
        Assertions.assertArrayEquals(new Object[] {"bigint", "YES"}, columnInfo("ddl_life", "quantity"));
        Assertions.assertEquals(10L, queryDirect("SELECT quantity FROM ddl_life WHERE id = 1").get(0)[0]);

        AlterTable setNotNull = alter("ddl_life", "setNullable", "name");
        setNotNull.nullable = Boolean.FALSE; // no NULLs present — allowed
        Assertions.assertNull(runMutation(setNotNull).errorMessage);
        Assertions.assertArrayEquals(new Object[] {"text", "NO"}, columnInfo("ddl_life", "name"));

        CreateIndex createIndex = new CreateIndex();
        createIndex.tableName = "ddl_life";
        createIndex.columns = Collections.singletonList("quantity");
        Assertions.assertNull(runMutation(createIndex).errorMessage);
        Assertions.assertTrue(indexExists("ddl_life", "ix_ddl_life_quantity")); // auto-named

        execDirect("CREATE TABLE ddl_child (id int PRIMARY KEY, parent_id int)");
        AddKey addFk = new AddKey();
        addFk.tableName = "ddl_child";
        addFk.keyType = "foreign";
        addFk.columns = Collections.singletonList("parent_id");
        addFk.refTable = "ddl_life";
        addFk.refColumns = Collections.singletonList("id");
        Assertions.assertNull(runMutation(addFk).errorMessage);
        Assertions.assertEquals(1, queryDirect("SELECT 1 FROM information_schema.table_constraints"
                + " WHERE table_name = 'ddl_child' AND constraint_name = 'fk_ddl_child_parent_id'"
                + " AND constraint_type = 'FOREIGN KEY'").size());

        DropIndex dropIndex = new DropIndex();
        dropIndex.tableName = "ddl_life";
        dropIndex.indexName = "ix_ddl_life_quantity";
        Assertions.assertNull(runMutation(dropIndex).errorMessage);
        Assertions.assertFalse(indexExists("ddl_life", "ix_ddl_life_quantity"));

        Assertions.assertNull(runMutation(dropTable("ddl_child", true)).errorMessage); // FK holder first
        Assertions.assertFalse(tableExists("ddl_child"));

        TruncateTable truncate = new TruncateTable();
        truncate.tableName = "ddl_life";
        truncate.confirmDestructive = true;
        result = runMutation(truncate);
        Assertions.assertNull(result.errorMessage);
        Assertions.assertEquals("truncate-table", result.plan.destructive.get(0).code);
        Assertions.assertEquals(2L, result.plan.destructive.get(0).count);
        Assertions.assertEquals(0, countDirect("ddl_life"));

        result = runMutation(dropTable("ddl_life", true));
        Assertions.assertNull(result.errorMessage);
        Assertions.assertEquals("drop-table", result.plan.destructive.get(0).code);
        Assertions.assertEquals(0L, result.plan.destructive.get(0).count); // truncated above
        Assertions.assertFalse(tableExists("ddl_life"));
    }

    @DisplayName("dryRun dropTable: schema unchanged (second connection), plan carries the exact SQL + live row count")
    @Test
    public void dryRun_dropTable_noChanges() throws Exception {
        execDirect("CREATE TABLE ddl_dry_drop (id int)");
        execDirect("INSERT INTO ddl_dry_drop VALUES (1), (2), (3)");
        DropTable m = dropTable("ddl_dry_drop", false);
        m.dryRun = true; // passes the same gates, executes nothing — no confirm needed to look
        MutationResult result = runMutation(m);
        Assertions.assertNull(result.errorMessage);
        Assertions.assertEquals(Collections.singletonList("DROP TABLE \"ddl_dry_drop\""), result.plan.statements);
        Assertions.assertEquals(Boolean.TRUE, result.plan.transactionalDdl);
        Assertions.assertEquals("drop-table", result.plan.destructive.get(0).code);
        Assertions.assertEquals("table", result.plan.destructive.get(0).path);
        Assertions.assertEquals(3L, result.plan.destructive.get(0).count);
        Assertions.assertTrue(tableExists("ddl_dry_drop"));
        Assertions.assertEquals(3, countDirect("ddl_dry_drop"));
    }

    @DisplayName("Destructive contract: dropTable without confirm refused with the recomputed plan, table intact; with confirm dropped")
    @Test
    public void dropTable_confirmContract() throws Exception {
        execDirect("CREATE TABLE ddl_confirm (id int)");
        execDirect("INSERT INTO ddl_confirm VALUES (1), (2)");
        DropTable refused = dropTable("ddl_confirm", false);
        MutationConfirmationRequiredException ex = Assertions.assertThrows(
                MutationConfirmationRequiredException.class, () -> runMutation(refused));
        Assertions.assertEquals("Postgres", ex.providerType);
        Assertions.assertEquals(Collections.singletonList("DROP TABLE \"ddl_confirm\""), ex.plan.statements);
        Assertions.assertEquals("drop-table", ex.plan.destructive.get(0).code);
        Assertions.assertEquals(2L, ex.plan.destructive.get(0).count);
        Assertions.assertTrue(tableExists("ddl_confirm"));
        Assertions.assertEquals(2, countDirect("ddl_confirm"));

        MutationResult result = runMutation(dropTable("ddl_confirm", true));
        Assertions.assertNull(result.errorMessage);
        Assertions.assertFalse(tableExists("ddl_confirm"));
    }

    @DisplayName("Hard refusal: setNullable(false) with NULLs present -> not-null-violation with the live count, no DDL ran")
    @Test
    public void setNullable_withNulls_hardRefusal() throws Exception {
        execDirect("CREATE TABLE ddl_nn (id int, name text)");
        execDirect("INSERT INTO ddl_nn VALUES (1, 'a'), (2, NULL)");
        AlterTable m = alter("ddl_nn", "setNullable", "name");
        m.nullable = Boolean.FALSE;
        m.confirmDestructive = true; // confirmation cannot help — the DDL would fail
        MutationValidationException ex = Assertions.assertThrows(MutationValidationException.class, () -> runMutation(m));
        Assertions.assertEquals("not-null-violation", ex.refusals.get(0).code);
        Assertions.assertEquals("columns.name", ex.refusals.get(0).path);
        Assertions.assertEquals(1L, ex.refusals.get(0).count);
        Assertions.assertArrayEquals(new Object[] {"text", "YES"}, columnInfo("ddl_nn", "name")); // untouched
    }

    @DisplayName("Hard refusal: NOT NULL addColumn without default refused on a populated table, allowed on an empty one")
    @Test
    public void addColumn_requiredWithoutDefault() throws Exception {
        execDirect("CREATE TABLE ddl_req (id int)");
        execDirect("INSERT INTO ddl_req VALUES (1)");
        AlterTable m = alter("ddl_req", "addColumn", null);
        m.column = spec("code", "string", false, null);
        MutationValidationException ex = Assertions.assertThrows(MutationValidationException.class, () -> runMutation(m));
        Assertions.assertEquals("required-without-default", ex.refusals.get(0).code);
        Assertions.assertEquals("columns.code", ex.refusals.get(0).path);
        Assertions.assertEquals(1L, ex.refusals.get(0).count);
        Assertions.assertNull(columnInfo("ddl_req", "code"));

        execDirect("CREATE TABLE ddl_req_empty (id int)");
        AlterTable allowed = alter("ddl_req_empty", "addColumn", null);
        allowed.column = spec("code", "string", false, null);
        Assertions.assertNull(runMutation(allowed).errorMessage); // pre-check count 0 — no refusal
        Assertions.assertArrayEquals(new Object[] {"text", "NO"}, columnInfo("ddl_req_empty", "code"));
    }

    @DisplayName("changeType over non-castable data: db-error reported in the payload, rolled back (PG transactional DDL)")
    @Test
    public void changeType_nonCastable_dbErrorRolledBack() throws Exception {
        execDirect("CREATE TABLE ddl_ct (id int, val text)");
        execDirect("INSERT INTO ddl_ct VALUES (1, 'abc')");
        AlterTable m = alter("ddl_ct", "changeType", "val");
        m.newType = "int";
        m.confirmDestructive = true;
        MutationResult result = runMutation(m);
        Assertions.assertNotNull(result.errorMessage);
        Assertions.assertEquals(1, result.errorCount);
        Assertions.assertArrayEquals(new Object[] {"text", "YES"}, columnInfo("ddl_ct", "val")); // type unchanged
        Assertions.assertEquals("abc", queryDirect("SELECT val FROM ddl_ct WHERE id = 1").get(0)[0]);
    }

    @DisplayName("Transaction honesty: CREATE TABLE + failing CREATE INDEX rolls the table back (PG transactional DDL)")
    @Test
    public void createTable_failingIndex_atomicRollback() throws Exception {
        CreateTable m = new CreateTable();
        m.tableName = "ddl_tx";
        m.columns = Collections.singletonList(spec("id", "int", true, null));
        m.indexes = Collections.singletonList(new IndexSpec());
        m.indexes.get(0).columns = Collections.singletonList("missing"); // valid identifier, no such column
        MutationResult result = runMutation(m);
        Assertions.assertNotNull(result.errorMessage);
        Assertions.assertFalse(tableExists("ddl_tx")); // the successful CREATE TABLE was rolled back too
    }

    @DisplayName("Capability gate: DDL against a supportsDdl=false provider is a structured capability refusal, before any connection")
    @Test
    public void ddl_onDdlLessProvider_capabilityRefused() {
        DropTable m = dropTable("ddl_life", true);
        m.connection = connection;
        FuncCall call = new FuncCall();
        call.func = m;
        call.options = new HashMap<>();
        Assertions.assertThrows(UnsupportedOperationException.class,
                () -> MutationRunner.execute(new ClickHouseProvider(), call));
    }

    @DisplayName("Single-op rule: MutationBatch with a DDL member or a mode=create insert is a structured validation refusal")
    @Test
    public void batch_withDdlOrCreateMode_refused() {
        MutationBatch withDdl = new MutationBatch();
        withDdl.operations = Arrays.asList(
                insertRows("ddl_life", Arrays.asList("id"), Arrays.asList("int"),
                        Collections.singletonList(Arrays.asList((Object) 1.0d))),
                dropTable("ddl_life", true));
        Assertions.assertThrows(MutationValidationException.class, () -> runMutation(withDdl));

        InsertRows createMode = insertRows("ddl_life", Arrays.asList("id"), Arrays.asList("int"), null);
        createMode.mode = "create";
        MutationBatch withCreate = new MutationBatch();
        withCreate.operations = Collections.singletonList(createMode);
        Assertions.assertThrows(MutationValidationException.class, () -> runMutation(withCreate));
    }

    @DisplayName("DML dryRun: exact statements returned, nothing executed (single op and batch)")
    @Test
    public void dmlDryRun_statementsOnly() throws Exception {
        execDirect("CREATE TABLE ddl_dml (id int, note text)");
        execDirect("INSERT INTO ddl_dml VALUES (1, 'seed')");
        InsertRows insert = insertRows("ddl_dml", Arrays.asList("id", "note"), Arrays.asList("int", "string"),
                Collections.singletonList(Arrays.asList((Object) 2.0d, "x")));
        insert.dryRun = true;
        MutationResult result = runMutation(insert);
        Assertions.assertNull(result.errorMessage);
        Assertions.assertEquals(Collections.singletonList("INSERT INTO \"ddl_dml\" (\"id\", \"note\") VALUES (?, ?)"),
                result.plan.statements);
        Assertions.assertNull(result.plan.transactionalDdl); // DML plan
        Assertions.assertEquals(1, countDirect("ddl_dml"));

        DeleteRows delete = new DeleteRows();
        delete.tableName = "ddl_dml";
        delete.allowFullTable = true;
        MutationBatch batch = new MutationBatch();
        batch.operations = Arrays.asList(
                insertRows("ddl_dml", Arrays.asList("id", "note"), Arrays.asList("int", "string"),
                        Collections.singletonList(Arrays.asList((Object) 3.0d, "y"))),
                delete);
        batch.dryRun = true;
        result = runMutation(batch);
        Assertions.assertNull(result.errorMessage);
        Assertions.assertEquals(2, result.plan.statements.size());
        Assertions.assertEquals("DELETE FROM \"ddl_dml\"", result.plan.statements.get(1));
        Assertions.assertEquals(1, countDirect("ddl_dml")); // nothing ran
    }
}

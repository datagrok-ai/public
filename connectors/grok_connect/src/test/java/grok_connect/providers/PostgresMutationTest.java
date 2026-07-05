package grok_connect.providers;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.FuncCallBuilder;
import grok_connect.providers.utils.Provider;
import grok_connect.table_mutation.DeleteRows;
import grok_connect.table_mutation.InsertRows;
import grok_connect.table_mutation.MutationBatch;
import grok_connect.table_mutation.MutationResult;
import grok_connect.table_mutation.MutationRunner;
import grok_connect.table_mutation.MutationValidationException;
import grok_connect.table_mutation.TableMutation;
import grok_connect.table_mutation.UpdateRows;
import grok_connect.table_mutation.UpsertRows;
import grok_connect.table_query.FieldPredicate;
import grok_connect.utils.PatternMatcher;
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
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TimeZone;
import java.util.stream.Collectors;

/**
 * Integration tests for GROK-20327: MutationRunner against a containerized Postgres —
 * insert/update/delete/batch happy paths, rollback on mid-batch error, empty-WHERE refusal,
 * dg-type coercion round trip (bigint beyond 2^53, datetime), and injection probes.
 */
class PostgresMutationTest extends ContainerizedProviderBaseTest {
    protected PostgresMutationTest() {
        super(Provider.POSTGRESQL);
    }

    /** Direct JDBC, bypassing grok_connect — verification must not depend on the code under test. */
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

    @AfterAll
    public void dropScratchTables() throws SQLException {
        for (String table : new String[] {"mut_ins", "mut_upd", "mut_del", "mut_batch", "mut_full", "mut_probe", "mut_inj", "mut_upsert"})
            execDirect("DROP TABLE IF EXISTS " + table);
    }

    private void runDdl(String query) throws Exception {
        FuncCall call = FuncCallBuilder.fromQuery(query);
        call.func.connection = connection;
        provider.execute(call);
    }

    private MutationResult runMutation(TableMutation mutation) throws Exception {
        mutation.connection = connection;
        FuncCall call = new FuncCall();
        call.func = mutation;
        call.options = new HashMap<>();
        return MutationRunner.execute(provider, call);
    }

    private FieldPredicate predicate(String field, String dataType, String pattern, String op, Object... values) {
        FieldPredicate clause = new FieldPredicate(field, pattern, dataType);
        Map<String, Object> options = new HashMap<>();
        options.put("expression", pattern);
        options.put("op", op);
        options.put("values", Arrays.stream(values).map(Object::toString).collect(Collectors.toList()));
        clause.matcher = new PatternMatcher(options, field);
        return clause;
    }

    private InsertRows insertRows(String table, List<String> columns, List<String> columnTypes, List<List<Object>> rows) {
        InsertRows m = new InsertRows();
        m.tableName = table;
        m.columns = columns;
        m.columnTypes = columnTypes;
        m.rows = rows;
        return m;
    }

    @DisplayName("Insert: dg-type coercion round trip (Double->int, string bigint beyond 2^53, datetime, nulls)")
    @Test
    public void insert_typeCoercionRoundTrip() throws Exception {
        runDdl("CREATE TABLE mut_ins (id int PRIMARY KEY, big_id bigint, val double precision, active boolean, note text, created timestamp)");
        InsertRows m = insertRows("mut_ins",
                Arrays.asList("id", "big_id", "val", "active", "note", "created"),
                Arrays.asList("int", "bigint", "double", "bool", "string", "datetime"),
                Arrays.asList(
                        // values shaped exactly as Gson delivers them: numbers as Double, big bigints as String
                        Arrays.asList((Object) 1.0d, "9007199254740995", 1.5d, true, "alpha", "2026-07-05T12:34:56.789Z"),
                        Arrays.asList((Object) 2.0d, "-9223372036854775808", null, false, null, null),
                        Arrays.asList((Object) 3.0d, 42.0d, -2.25d, true, "gamma", "2026-01-01T00:00:00.000Z")));
        MutationResult result = runMutation(m);
        Assertions.assertNull(result.errorMessage);
        Assertions.assertEquals(3, result.affectedRows);

        try (Connection c = DriverManager.getConnection(container.getJdbcUrl(), container.getUsername(), container.getPassword());
             PreparedStatement statement = c.prepareStatement("SELECT big_id, val, active, note, created FROM mut_ins ORDER BY id");
             ResultSet rs = statement.executeQuery()) {
            Calendar utc = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
            Assertions.assertTrue(rs.next());
            Assertions.assertEquals(9007199254740995L, rs.getLong("big_id")); // beyond 2^53 — exact
            Assertions.assertEquals(1.5d, rs.getDouble("val"));
            Assertions.assertTrue(rs.getBoolean("active"));
            Assertions.assertEquals("alpha", rs.getString("note"));
            Assertions.assertEquals(Instant.parse("2026-07-05T12:34:56.789Z").toEpochMilli(),
                    rs.getTimestamp("created", utc).getTime());
            Assertions.assertTrue(rs.next());
            Assertions.assertEquals(Long.MIN_VALUE, rs.getLong("big_id"));
            rs.getDouble("val");
            Assertions.assertTrue(rs.wasNull());
            Assertions.assertNull(rs.getString("note"));
            Assertions.assertNull(rs.getTimestamp("created", utc));
            Assertions.assertTrue(rs.next());
            Assertions.assertEquals(42L, rs.getLong("big_id")); // small bigint may travel as a JSON number
            Assertions.assertEquals(Instant.parse("2026-01-01T00:00:00.000Z").toEpochMilli(),
                    rs.getTimestamp("created", utc).getTime());
            Assertions.assertFalse(rs.next());
        }
    }

    @DisplayName("Update: string 'contains' + numeric range predicates hit only the matching row")
    @Test
    public void update_withPatternPredicates() throws Exception {
        runDdl("CREATE TABLE mut_upd (id int PRIMARY KEY, name text, qty int, status text)");
        execDirect("INSERT INTO mut_upd VALUES (1, 'apple', 5, 'new'), (2, 'grape', 7, 'new'), (3, 'apple pie', 20, 'new')");
        UpdateRows m = new UpdateRows();
        m.tableName = "mut_upd";
        m.setColumns = Arrays.asList("status");
        m.setValues = Arrays.asList((Object) "matched");
        m.setTypes = Arrays.asList("string");
        m.whereClauses = new ArrayList<>();
        m.whereClauses.add(predicate("name", "string", "contains apple", "contains", "apple"));
        m.whereClauses.add(predicate("qty", "int", "1-10", "-", 1, 10));
        MutationResult result = runMutation(m);
        Assertions.assertNull(result.errorMessage);
        Assertions.assertEquals(1, result.affectedRows);
        Assertions.assertEquals("matched", queryDirect("SELECT status FROM mut_upd WHERE id = 1").get(0)[0]);
        Assertions.assertEquals(2, ((Number) queryDirect("SELECT count(*) FROM mut_upd WHERE status = 'new'").get(0)[0]).intValue());
    }

    @DisplayName("Delete: by predicate removes only matching rows")
    @Test
    public void delete_byPredicate() throws Exception {
        runDdl("CREATE TABLE mut_del (id int PRIMARY KEY, status text)");
        execDirect("INSERT INTO mut_del VALUES (1, 'obsolete'), (2, 'active'), (3, 'obsolete')");
        DeleteRows m = new DeleteRows();
        m.tableName = "mut_del";
        m.whereClauses = new ArrayList<>();
        m.whereClauses.add(predicate("status", "string", "obsolete", "equals", "obsolete"));
        MutationResult result = runMutation(m);
        Assertions.assertNull(result.errorMessage);
        Assertions.assertEquals(2, result.affectedRows);
        Assertions.assertEquals(1, countDirect("mut_del"));
        Assertions.assertEquals(2, queryDirect("SELECT id FROM mut_del").get(0)[0]);
    }

    @DisplayName("Batch: happy path executes all ops in order with perStatement counts")
    @Test
    public void batch_happyPath() throws Exception {
        runDdl("CREATE TABLE IF NOT EXISTS mut_batch (id int PRIMARY KEY, qty int)");
        execDirect("DELETE FROM mut_batch");
        MutationBatch batch = new MutationBatch();
        batch.tableName = "mut_batch";
        InsertRows insert = insertRows("mut_batch",
                Arrays.asList("id", "qty"), Arrays.asList("int", "int"),
                Arrays.asList(Arrays.asList((Object) 1.0d, 10.0d), Arrays.asList((Object) 2.0d, 20.0d)));
        UpdateRows update = new UpdateRows();
        update.tableName = "mut_batch";
        update.setColumns = Arrays.asList("qty");
        update.setValues = Arrays.asList((Object) 99.0d);
        update.setTypes = Arrays.asList("int");
        update.whereClauses = new ArrayList<>();
        update.whereClauses.add(predicate("id", "int", "1", "=", 1));
        DeleteRows delete = new DeleteRows();
        delete.tableName = "mut_batch";
        delete.whereClauses = new ArrayList<>();
        delete.whereClauses.add(predicate("id", "int", "2", "=", 2));
        batch.operations = Arrays.asList(insert, update, delete);
        MutationResult result = runMutation(batch);
        Assertions.assertNull(result.errorMessage);
        Assertions.assertEquals(4, result.affectedRows); // 2 inserted + 1 updated + 1 deleted
        Assertions.assertEquals(3, result.perStatement.size());
        Assertions.assertEquals(2, result.perStatement.get(0).affectedRows);
        Assertions.assertEquals(1, result.perStatement.get(1).affectedRows);
        Assertions.assertEquals(1, result.perStatement.get(2).affectedRows);
        Assertions.assertEquals(99, queryDirect("SELECT qty FROM mut_batch WHERE id = 1").get(0)[0]);
        Assertions.assertEquals(1, countDirect("mut_batch"));
        execDirect("DELETE FROM mut_batch");
    }

    @DisplayName("Batch: mid-batch constraint violation rolls back all ops and names the failing statement")
    @Test
    public void batch_rollbackOnMidBatchError() throws Exception {
        runDdl("CREATE TABLE IF NOT EXISTS mut_batch (id int PRIMARY KEY, qty int)");
        execDirect("DELETE FROM mut_batch");
        MutationBatch batch = new MutationBatch();
        batch.tableName = "mut_batch";
        InsertRows first = insertRows("mut_batch",
                Arrays.asList("id", "qty"), Arrays.asList("int", "int"),
                Arrays.asList(Arrays.asList((Object) 10.0d, 1.0d)));
        UpdateRows second = new UpdateRows();
        second.tableName = "mut_batch";
        second.setColumns = Arrays.asList("qty");
        second.setValues = Arrays.asList((Object) 2.0d);
        second.setTypes = Arrays.asList("int");
        second.whereClauses = new ArrayList<>();
        second.whereClauses.add(predicate("id", "int", "10", "=", 10));
        InsertRows third = insertRows("mut_batch", // PK violation
                Arrays.asList("id", "qty"), Arrays.asList("int", "int"),
                Arrays.asList(Arrays.asList((Object) 10.0d, 3.0d)));
        batch.operations = Arrays.asList(first, second, third);
        MutationResult result = runMutation(batch);
        Assertions.assertNotNull(result.errorMessage);
        Assertions.assertEquals(1, result.errorCount);
        Assertions.assertEquals(2, result.errors.get(0).index); // failing statement index
        Assertions.assertEquals(0, countDirect("mut_batch")); // ops 1-2 rolled back
    }

    @DisplayName("Delete: empty WHERE refused without allowFullTable, wipes the table with it")
    @Test
    public void delete_emptyWhere() throws Exception {
        runDdl("CREATE TABLE mut_full (id int PRIMARY KEY)");
        execDirect("INSERT INTO mut_full VALUES (1), (2), (3)");
        DeleteRows refused = new DeleteRows();
        refused.tableName = "mut_full";
        Assertions.assertThrows(MutationValidationException.class, () -> runMutation(refused));
        Assertions.assertEquals(3, countDirect("mut_full"));
        DeleteRows allowed = new DeleteRows();
        allowed.tableName = "mut_full";
        allowed.allowFullTable = true;
        MutationResult result = runMutation(allowed);
        Assertions.assertNull(result.errorMessage);
        Assertions.assertEquals(3, result.affectedRows);
        Assertions.assertEquals(0, countDirect("mut_full"));
    }

    @DisplayName("Injection probes: hostile value stays a bound literal, hostile column name cannot escape brackets")
    @Test
    public void injection_probes() throws Exception {
        runDdl("CREATE TABLE mut_probe (id int)");
        runDdl("CREATE TABLE mut_inj (id int, note text)");
        String hostileValue = "'); DROP TABLE mut_probe; --";
        MutationResult valueProbe = runMutation(insertRows("mut_inj",
                Arrays.asList("id", "note"), Arrays.asList("int", "string"),
                Arrays.asList(Arrays.asList((Object) 1.0d, hostileValue))));
        Assertions.assertNull(valueProbe.errorMessage);
        Assertions.assertEquals(hostileValue, queryDirect("SELECT note FROM mut_inj WHERE id = 1").get(0)[0]);
        Assertions.assertEquals(0, countDirect("mut_probe")); // still exists

        // hostile column identifier is refused at the mutation boundary (structured validation,
        // not a downstream db-error): SQL is never built, nothing executes
        Assertions.assertThrows(MutationValidationException.class, () -> runMutation(insertRows("mut_inj",
                Arrays.asList("id", "note\"; DROP TABLE mut_probe; --"), Arrays.asList("int", "string"),
                Arrays.asList(Arrays.asList((Object) 2.0d, "x")))));
        Assertions.assertEquals(0, countDirect("mut_probe")); // still exists
        Assertions.assertEquals(1, countDirect("mut_inj")); // probe insert never happened
    }

    @DisplayName("Upsert: seed 3, merge a 4-row payload (2 update + 2 insert by matchKeys), final state correct")
    @Test
    public void upsert_seedAndMerge() throws Exception {
        runDdl("CREATE TABLE mut_upsert (id int PRIMARY KEY, region text, amount double precision)");
        execDirect("INSERT INTO mut_upsert VALUES (1, 'east', 10), (2, 'west', 20), (3, 'north', 30)");
        UpsertRows m = new UpsertRows();
        m.tableName = "mut_upsert";
        m.columns = Arrays.asList("id", "region", "amount");
        m.columnTypes = Arrays.asList("int", "string", "double");
        m.matchKeys = Arrays.asList("id");
        m.rows = Arrays.asList(
                Arrays.asList((Object) 1.0d, "east", 111.0d),   // update
                Arrays.asList((Object) 2.0d, "west", 222.0d),   // update
                Arrays.asList((Object) 4.0d, "south", 40.0d),   // insert
                Arrays.asList((Object) 5.0d, "central", 50.0d));// insert
        MutationResult result = runMutation(m);
        Assertions.assertNull(result.errorMessage);
        Assertions.assertEquals(4, result.affectedRows); // Postgres reports 1 per affected row
        Assertions.assertEquals(5, countDirect("mut_upsert"));
        Assertions.assertEquals(111.0d, ((Number) queryDirect("SELECT amount FROM mut_upsert WHERE id = 1").get(0)[0]).doubleValue());
        Assertions.assertEquals(222.0d, ((Number) queryDirect("SELECT amount FROM mut_upsert WHERE id = 2").get(0)[0]).doubleValue());
        Assertions.assertEquals(30.0d, ((Number) queryDirect("SELECT amount FROM mut_upsert WHERE id = 3").get(0)[0]).doubleValue()); // untouched
        Assertions.assertEquals("south", queryDirect("SELECT region FROM mut_upsert WHERE id = 4").get(0)[0]);
        Assertions.assertEquals("central", queryDirect("SELECT region FROM mut_upsert WHERE id = 5").get(0)[0]);
    }

    @DisplayName("Upsert: absent matchKeys is a structured validation error")
    @Test
    public void upsert_noMatchKeys_refused() {
        UpsertRows m = new UpsertRows();
        m.tableName = "mut_upsert";
        m.columns = Arrays.asList("id");
        m.columnTypes = Arrays.asList("int");
        m.rows = Arrays.asList(Arrays.asList((Object) 1.0d));
        Assertions.assertThrows(grok_connect.table_mutation.MutationValidationException.class, () -> runMutation(m));
    }
}

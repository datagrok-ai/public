package grok_connect.table_mutation;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataQueryRunResult;
import grok_connect.connectors_info.FuncCall;
import grok_connect.connectors_info.FuncParam;
import grok_connect.providers.MsSqlDataProvider;
import grok_connect.providers.PostgresDataProvider;
import grok_connect.table_query.FieldPredicate;
import grok_connect.table_query.TableQuery;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.PatternMatcher;
import grok_connect.utils.ProviderManager;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInstance;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * SQL-emission tests for table_mutation (GROK-20327) — no Docker. Pins the exact SQL for
 * insert/update/delete against Postgres ("") and the generic descriptor brackets ([]),
 * the empty-WHERE delete safety rail, predicate-compilation parity with TableQuery.toSql,
 * and the mutations-never-enter-the-query-path guards.
 */
@TestInstance(TestInstance.Lifecycle.PER_CLASS)
public class MutationSqlTest {
    private final PostgresDataProvider postgres = new PostgresDataProvider();
    private final MsSqlDataProvider generic = new MsSqlDataProvider(); // default "[]" nameBrackets

    @BeforeAll
    public void init() {
        if (GrokConnect.providerManager == null)
            GrokConnect.providerManager = new ProviderManager();
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

    private InsertRows insertOrders() {
        InsertRows m = new InsertRows();
        m.tableName = "orders";
        m.schema = "public";
        m.columns = Arrays.asList("id", "note");
        m.columnTypes = Arrays.asList("int", "string");
        return m;
    }

    private String lines(String... lines) {
        return String.join(System.lineSeparator(), lines);
    }

    @DisplayName("insertSql: Postgres quoting, one row of placeholders")
    @Test
    public void insertSql_postgres() {
        Assertions.assertEquals("INSERT INTO \"public\".\"orders\" (\"id\", \"note\") VALUES (?, ?)",
                postgres.insertSql(insertOrders()));
    }

    @DisplayName("insertSql: generic [] descriptor brackets")
    @Test
    public void insertSql_genericBrackets() {
        Assertions.assertEquals("INSERT INTO [public].[orders] ([id], [note]) VALUES (?, ?)",
                generic.insertSql(insertOrders()));
    }

    @DisplayName("insertSql: empty columns list refused")
    @Test
    public void insertSql_noColumns_refused() {
        InsertRows m = insertOrders();
        m.columns = new ArrayList<>();
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.insertSql(m));
    }

    @DisplayName("updateSql: SET placeholders + predicate WHERE, params collected")
    @Test
    public void updateSql_postgres() {
        UpdateRows m = new UpdateRows();
        m.tableName = "orders";
        m.schema = "public";
        m.setColumns = Arrays.asList("status", "qty");
        m.setValues = Arrays.asList("shipped", 5);
        m.setTypes = Arrays.asList("string", "int");
        m.whereClauses = new ArrayList<>();
        m.whereClauses.add(predicate("id", "int", ">3", ">", 3));
        List<FuncParam> params = new ArrayList<>();
        Assertions.assertEquals(lines(
                "UPDATE \"public\".\"orders\" SET \"status\" = ?, \"qty\" = ?",
                "WHERE",
                "  ((\"id\" > @id))"),
                postgres.updateSql(m, params));
        Assertions.assertEquals(1, params.size());
        Assertions.assertEquals("id", params.get(0).name);
        Assertions.assertEquals("int", params.get(0).propertyType);
        Assertions.assertEquals(3, params.get(0).value);
    }

    @DisplayName("updateSql: parallel-list size mismatch refused")
    @Test
    public void updateSql_parallelListMismatch_refused() {
        UpdateRows m = new UpdateRows();
        m.tableName = "orders";
        m.setColumns = Arrays.asList("status", "qty");
        m.setValues = Arrays.asList((Object) "shipped");
        m.setTypes = Arrays.asList("string", "int");
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.updateSql(m, new ArrayList<>()));
    }

    @DisplayName("deleteSql: predicate WHERE with pattern converter")
    @Test
    public void deleteSql_postgres() {
        DeleteRows m = new DeleteRows();
        m.tableName = "orders";
        m.schema = "public";
        m.whereClauses = new ArrayList<>();
        m.whereClauses.add(predicate("status", "string", "contains obsolete", "contains", "obsolete"));
        List<FuncParam> params = new ArrayList<>();
        Assertions.assertEquals(lines(
                "DELETE FROM \"public\".\"orders\"",
                "WHERE",
                "  ((LOWER(\"status\") LIKE @status))"),
                postgres.deleteSql(m, params));
        Assertions.assertEquals(1, params.size());
        Assertions.assertEquals("%obsolete%", params.get(0).value);
    }

    @DisplayName("deleteSql: empty WHERE without allowFullTable is refused; with the flag it emits")
    @Test
    public void deleteSql_emptyWhere() {
        DeleteRows m = new DeleteRows();
        m.tableName = "orders";
        m.schema = "public";
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.deleteSql(m, new ArrayList<>()));
        m.allowFullTable = true;
        Assertions.assertEquals("DELETE FROM \"public\".\"orders\"", postgres.deleteSql(m, new ArrayList<>()));
    }

    @DisplayName("deleteSql: whereOp is restricted to and/or")
    @Test
    public void deleteSql_whereOpInjection_refused() {
        DeleteRows m = new DeleteRows();
        m.tableName = "orders";
        m.whereOp = "and (1=1); DROP TABLE orders; --";
        m.whereClauses = new ArrayList<>();
        m.whereClauses.add(predicate("id", "int", ">3", ">", 3));
        m.whereClauses.add(predicate("qty", "int", ">1", ">", 1));
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.deleteSql(m, new ArrayList<>()));
    }

    @DisplayName("Predicate parity: identical WHERE through TableQuery.toSql and deleteSql")
    @Test
    public void predicateParity_tableQueryAndDelete() {
        TableQuery query = new TableQuery();
        query.tableName = "orders";
        query.schema = "public";
        query.connection = new DataConnection();
        query.connection.dataSource = postgres.descriptor.type;
        query.params = new ArrayList<>();
        query.whereClauses.add(predicate("qty", "int", "5-10", "-", 5, 10));
        query.whereClauses.add(predicate("status", "string", "contains new", "contains", "new"));
        String querySql = postgres.queryTableSql(query);

        DeleteRows delete = new DeleteRows();
        delete.tableName = "orders";
        delete.schema = "public";
        delete.whereClauses = new ArrayList<>();
        delete.whereClauses.add(predicate("qty", "int", "5-10", "-", 5, 10));
        delete.whereClauses.add(predicate("status", "string", "contains new", "contains", "new"));
        List<FuncParam> params = new ArrayList<>();
        String deleteSql = postgres.deleteSql(delete, params);

        String queryWhere = querySql.substring(querySql.indexOf("WHERE"));
        String deleteWhere = deleteSql.substring(deleteSql.indexOf("WHERE"));
        Assertions.assertEquals(queryWhere, deleteWhere);
        Assertions.assertEquals(3, params.size()); // qtyR0, qtyR1, status
    }

    @DisplayName("Identifier validation: hostile INSERT column refused before SQL is built")
    @Test
    public void insertSql_hostileColumn_refused() {
        InsertRows m = insertOrders();
        m.columns = Arrays.asList("id", "note\"; DROP TABLE orders; --");
        m.columnTypes = Arrays.asList("int", "string");
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.insertSql(m));
    }

    @DisplayName("Identifier validation: hostile UPDATE setColumn refused")
    @Test
    public void updateSql_hostileSetColumn_refused() {
        UpdateRows m = new UpdateRows();
        m.tableName = "orders";
        m.setColumns = Arrays.asList("qty\" = 0; DROP TABLE orders; --");
        m.setValues = Arrays.asList((Object) 1);
        m.setTypes = Arrays.asList("int");
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.updateSql(m, new ArrayList<>()));
    }

    @DisplayName("Identifier validation: hostile predicate field refused")
    @Test
    public void deleteSql_hostilePredicateField_refused() {
        DeleteRows m = new DeleteRows();
        m.tableName = "orders";
        m.whereClauses = new ArrayList<>();
        m.whereClauses.add(predicate("id\") OR 1=1 --", "int", ">0", ">", 0));
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.deleteSql(m, new ArrayList<>()));
    }

    @DisplayName("Identifier validation: leading-quote identifier refused (addBrackets would pass it verbatim)")
    @Test
    public void insertSql_leadingQuoteIdentifier_refused() {
        InsertRows m = insertOrders();
        m.columns = Arrays.asList("id", "\"note");
        m.columnTypes = Arrays.asList("int", "string");
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.insertSql(m));
        // and a hostile table name is refused too
        InsertRows table = insertOrders();
        table.tableName = "orders\"; DROP TABLE orders; --";
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.insertSql(table));
    }

    @DisplayName("Identifier validation: dotted schema.table and table.column stay valid")
    @Test
    public void validateIdentifier_dottedNamesValid() {
        InsertRows m = insertOrders();
        m.tableName = "public.orders";
        Assertions.assertEquals("INSERT INTO \"public\".\"orders\" (\"id\", \"note\") VALUES (?, ?)",
                postgres.insertSql(m));
        DeleteRows d = new DeleteRows();
        d.tableName = "orders";
        d.whereClauses = new ArrayList<>();
        d.whereClauses.add(predicate("orders.id", "int", ">0", ">", 0));
        Assertions.assertDoesNotThrow(() -> postgres.deleteSql(d, new ArrayList<>()));
    }

    @DisplayName("QueryManager refuses mutation FuncCalls (dryRun/streaming path exclusion)")
    @Test
    public void queryManager_refusesMutations() {
        InsertRows m = insertOrders();
        m.connection = new DataConnection();
        m.connection.dataSource = postgres.descriptor.type;
        FuncCall call = new FuncCall();
        call.func = m;
        call.options = new HashMap<>();
        String json = GrokConnect.gson.toJson(call);
        MutationValidationException e = Assertions.assertThrows(MutationValidationException.class,
                () -> new grok_connect.utils.QueryManager(json));
        Assertions.assertTrue(e.getMessage().contains("/mutate"));
    }

    @DisplayName("JdbcDataProvider.execute refuses mutation FuncCalls (POST /query path exclusion)")
    @Test
    public void execute_refusesMutations() {
        DeleteRows m = new DeleteRows();
        m.tableName = "orders";
        m.allowFullTable = true;
        FuncCall call = new FuncCall();
        call.func = m;
        call.options = new HashMap<>();
        GrokConnectException e = Assertions.assertThrows(GrokConnectException.class, () -> postgres.execute(call));
        Assertions.assertTrue(e.getMessage().contains("/mutate"));
    }

    @DisplayName("POST /query mutation refusal surfaces the structured guard message, not an opaque 500")
    @Test
    public void queryPath_mutationRefusal_structuredMessage() {
        DeleteRows m = new DeleteRows();
        m.tableName = "orders";
        m.allowFullTable = true;
        FuncCall call = new FuncCall();
        call.func = m;
        call.options = new HashMap<>();
        GrokConnectException ex = Assertions.assertThrows(GrokConnectException.class, () -> postgres.execute(call));
        // reproduce the fixed POST /query catch: a GrokConnectException with a null cause must
        // pack itself (not its null cause, which previously NPE'd printError -> generic 500)
        DataQueryRunResult result = new DataQueryRunResult();
        Exception forPack = ex.getClass().equals(GrokConnectException.class) && ex.getCause() != null
                ? (Exception) ex.getCause() : ex;
        Assertions.assertDoesNotThrow(() -> GrokConnect.packException(result, forPack));
        Assertions.assertNotNull(result.errorMessage);
        Assertions.assertTrue(result.errorMessage.contains("/mutate"));
        // latent NPE closed: printError tolerates a null throwable
        Assertions.assertDoesNotThrow(() -> GrokConnect.printError(null));
    }
}

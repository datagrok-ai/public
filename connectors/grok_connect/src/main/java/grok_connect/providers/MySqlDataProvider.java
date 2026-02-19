package grok_connect.providers;

import java.sql.*;
import java.util.*;
import java.util.stream.Collectors;
import grok_connect.managers.ColumnManager;
import grok_connect.managers.bool_column.MySqlMssqlBoolColumnManager;
import grok_connect.connectors_info.*;
import grok_connect.resultset.DefaultResultSetManager;
import grok_connect.resultset.ResultSetManager;
import grok_connect.table_query.AggrFunctionInfo;
import grok_connect.table_query.Stats;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.QueryCancelledByUser;
import serialization.DataFrame;
import serialization.StringColumn;
import serialization.Types;

public class MySqlDataProvider extends JdbcDataProvider {

    public MySqlDataProvider() {
        driverClassName = "com.mysql.cj.jdbc.Driver";
        descriptor = new DataSource();
        descriptor.type = "MySQL";
        descriptor.description = "Query MySQL database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(DbCredentials.getSsl());
        descriptor.credentialsTemplate = DbCredentials.getDbCredentialsTemplate();
        descriptor.canBrowseSchema = true;
        descriptor.nameBrackets = "`";
        descriptor.commentStart = "-- ";

        descriptor.typesMap = new HashMap<String, String>() {{
            put("bool", Types.BOOL);
            put("boolean", Types.BOOL);
            put("#bit(1)", Types.BOOL);
            put("#bit.*", Types.BIG_INT);
            put("#.*int", Types.INT);
            put("bigint", Types.BIG_INT);
            put("decimal", Types.FLOAT);
            put("float", Types.FLOAT);
            put("double", Types.FLOAT);
            put("double precision", Types.FLOAT);
            put("char", Types.STRING);
            put("varchar", Types.STRING);
            put("#.*text", Types.STRING);
            put("date", Types.DATE_TIME);
            put("datetime", Types.DATE_TIME);
            put("timestamp", Types.DATE_TIME);
            put("time", Types.DATE_TIME);
            put("year", Types.DATE_TIME);
            put("binary", Types.BLOB);
            put("varbinary", Types.BLOB);
            put("geometry", Types.OBJECT);
            put("point", Types.OBJECT);
            put("json", Types.OBJECT);
        }};
        descriptor.aggregations.add(new AggrFunctionInfo(Stats.STDEV, "std(#)", Types.dataFrameNumericTypes));
    }

    @Override
    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString()) {
            properties.setProperty("zeroDateTimeBehavior", "convertToNull");
            if (conn.ssl()) {
                properties.setProperty("useSSL", "true");
                properties.setProperty("verifyServerCertificate", "false");
            }
        }
        properties.setProperty("socketTimeout", "180000");
        return properties;
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:mysql://" + conn.getServer() + port + "/" + conn.getDb();
    }

    @Override
    public DataFrame getSchemas(DataConnection connection) throws QueryCancelledByUser, GrokConnectException {
        FuncCall queryRun = new FuncCall();
        queryRun.func = new DataQuery();
        queryRun.func.query = "SELECT schema_name AS table_schema FROM information_schema.schemata ORDER BY schema_name";
        queryRun.func.connection = connection;
        return execute(queryRun);
    }
    @Override
    public String getSchemaSql(String db, String schema, String table, boolean includeKeyInfo) {
        List<String> filters = new ArrayList<>();

        if (db != null && db.length() != 0)
            filters.add("c.table_schema = '" + db + "'");

        if (table != null)
            filters.add("(c.table_name = '" + table + "')");

        String whereClause = filters.size() != 0 ? "WHERE " + String.join(" AND \n", filters) : "";

        StringBuilder sql = new StringBuilder();

        sql.append("SELECT ")
                .append("c.table_schema AS table_schema, ")
                .append("c.table_name AS table_name, ")
                .append("c.column_name AS column_name, ")
                .append("c.data_type AS data_type, ")
                .append("CASE t.table_type WHEN 'VIEW' THEN 1 ELSE 0 END AS is_view");

        if (includeKeyInfo) {
            sql.append(", ")
                    .append("CASE WHEN pk.column_name IS NOT NULL THEN 1 ELSE 0 END AS is_primary_key, ")
                    .append("CASE WHEN pk.column_name IS NOT NULL OR uq.column_name IS NOT NULL THEN 1 ELSE 0 END AS is_unique ");
        }

        sql.append(" FROM information_schema.columns c ")
                .append("JOIN information_schema.tables t ")
                .append("  ON t.table_name = c.table_name ")
                .append(" AND t.table_schema = c.table_schema ")
                .append(" AND t.table_catalog = c.table_catalog ");

        if (includeKeyInfo) {
            sql.append("LEFT JOIN ( ")
                    .append("    SELECT kcu.table_schema, kcu.table_name, kcu.column_name ")
                    .append("    FROM information_schema.table_constraints tc ")
                    .append("    JOIN information_schema.key_column_usage kcu ")
                    .append("      ON tc.constraint_name = kcu.constraint_name ")
                    .append("     AND tc.table_schema = kcu.table_schema ")
                    .append("    WHERE tc.constraint_type = 'PRIMARY KEY' ")
                    .append(") pk ON pk.table_schema = c.table_schema ")
                    .append("     AND pk.table_name = c.table_name ")
                    .append("     AND pk.column_name = c.column_name ");

            sql.append("LEFT JOIN ( ")
                    .append("    SELECT kcu.table_schema, kcu.table_name, kcu.column_name ")
                    .append("    FROM information_schema.table_constraints tc ")
                    .append("    JOIN information_schema.key_column_usage kcu ")
                    .append("      ON tc.constraint_name = kcu.constraint_name ")
                    .append("     AND tc.table_schema = kcu.table_schema ")
                    .append("    WHERE tc.constraint_type = 'UNIQUE' ")
                    .append(") uq ON uq.table_schema = c.table_schema ")
                    .append("     AND uq.table_name = c.table_name ")
                    .append("     AND uq.column_name = c.column_name ");
        }

        sql.append(" ").append(whereClause);
        sql.append(" ORDER BY c.ORDINAL_POSITION;");
        return sql.toString();
    }


    @Override
    protected String getRegexQuery(String columnName, String regexExpression) {
        return String.format("%s REGEXP '%s'", columnName, regexExpression);
    }

    @Override
    protected void appendQueryParam(DataQuery dataQuery, String paramName, StringBuilder queryBuffer) {
        FuncParam param = dataQuery.getParam(paramName);
        if (param.propertyType.equals("list")) {
            @SuppressWarnings("unchecked")
            List<String> values = ((ArrayList<String>) param.value);
            queryBuffer.append(values.stream().map(value -> "?").collect(Collectors.joining(", ")));
        } else {
            queryBuffer.append("?");
        }
    }

    @Override
    protected int setArrayParamValue(PreparedStatement statement, int n, FuncParam param) throws SQLException {
        @SuppressWarnings (value="unchecked")
        ArrayList<Object> lst = (ArrayList<Object>)param.value;
        if (lst == null || lst.size() == 0) {
            statement.setObject(n, null);
            return 0;
        }
        for (int i = 0; i < lst.size(); i++) {
            statement.setObject(n + i, lst.get(i));
        }
        return lst.size() - 1;
    }

    @Override
    public ResultSetManager getResultSetManager() {
        Map<String, ColumnManager<?>> defaultManagersMap = DefaultResultSetManager.getDefaultManagersMap();
        defaultManagersMap.put(Types.BOOL, new MySqlMssqlBoolColumnManager());
        return DefaultResultSetManager.fromManagersMap(defaultManagersMap);
    }

    @Override
    public DataFrame getForeignKeys(DataConnection conn, String schema) throws GrokConnectException {
        try (Connection connection = getConnection(conn)) {
            DatabaseMetaData meta = connection.getMetaData();
            List<String> tables = new ArrayList<>();
            try (ResultSet tablesRs = meta.getTables(conn.getDb(), null, null, new String[]{"TABLE", "VIEW"})) {
                while (tablesRs.next())
                    tables.add(tablesRs.getString("TABLE_NAME"));
            }

            DataFrame result = DataFrame.fromColumns(new StringColumn("table_schema"),
                    new StringColumn("constraint_name"), new StringColumn("table_name"),
                    new StringColumn("column_name"), new StringColumn("foreign_table_name"), new StringColumn("foreign_column_name"));
            if (!tables.isEmpty()) {
                for (String t : tables)
                    try (ResultSet info = meta.getExportedKeys(conn.getDb(), null, t)) {
                        while(info.next())
                            result.addRow(info.getString("FKTABLE_SCHEM"), info.getString("FK_NAME"),
                                    info.getString("FKTABLE_NAME"), info.getString("FKCOLUMN_NAME"),
                                    info.getString("PKTABLE_NAME"), info.getString("PKCOLUMN_NAME"));
                    }
            }
            return result;
        } catch (SQLException e) {
            throw new GrokConnectException(e);
        }
    }

    @Override
    public String getCommentsQuery(DataConnection connection) {
        return "--input: string schema\n" +
                "(\n" +
                "    SELECT\n" +
                "        t.table_schema AS table_schema,\n" +
                "        CASE t.table_type\n" +
                "            WHEN 'VIEW' THEN 'view'\n" +
                "            ELSE 'table'\n" +
                "        END AS object_type,\n" +
                "        t.table_name AS table_name,\n" +
                "        NULL AS column_name,\n" +
                "        t.table_comment AS comment\n" +
                "    FROM information_schema.tables t\n" +
                "    WHERE t.table_schema = @schema\n" +
                ")\n" +
                "\n" +
                "UNION ALL\n" +
                "\n" +
                "-- COLUMN comments\n" +
                "(\n" +
                "    SELECT\n" +
                "        c.table_schema AS table_schema,\n" +
                "        'column' AS object_type,\n" +
                "        c.table_name AS table_name,\n" +
                "        c.column_name AS column_name,\n" +
                "        c.column_comment AS comment\n" +
                "    FROM information_schema.columns c\n" +
                "    WHERE c.table_schema = @schema\n" +
                ")\n" +
                "\n" +
                "ORDER BY object_type, table_name, column_name;\n";
    }
}

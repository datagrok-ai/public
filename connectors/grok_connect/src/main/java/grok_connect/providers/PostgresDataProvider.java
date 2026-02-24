package grok_connect.providers;

import java.sql.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Properties;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncParam;
import grok_connect.table_query.AggrFunctionInfo;
import grok_connect.table_query.Stats;
import grok_connect.utils.PatternMatcher;
import grok_connect.utils.Prop;
import grok_connect.utils.Property;
import org.postgresql.util.PGobject;
import serialization.Types;

public class PostgresDataProvider extends JdbcDataProvider {
    public PostgresDataProvider() {
        driverClassName = "org.postgresql.Driver";

        descriptor = new DataSource();
        descriptor.type = "Postgres";
        descriptor.description = "Query Postgres database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(DbCredentials.getSsl());
        descriptor.credentialsTemplate = DbCredentials.getDbCredentialsTemplate();
        descriptor.nameBrackets = "\"";

        descriptor.canBrowseSchema = true;
        descriptor.supportCatalogs = true;
        descriptor.defaultSchema = "public";
        descriptor.typesMap = new HashMap<String, String>() {{
            put("smallint", Types.INT);
            put("int", Types.INT);
            put("integer", Types.INT);
            put("bigint", Types.BIG_INT);
            put("real", Types.FLOAT);
            put("double precision", Types.FLOAT);
            put("numeric", Types.FLOAT);
            put("#character.*", Types.STRING);
            put("#varchar.*", Types.STRING);
            put("text", Types.STRING);
            put("boolean", Types.BOOL);
            put("date", Types.DATE_TIME);
            put("#timestamp.*", Types.DATE_TIME);
            put("cidr", Types.STRING);
            put("ARRAY", Types.LIST);
            put("USER_DEFINED", Types.STRING);
            put("bit.*", Types.BIG_INT);
            put("uuid", Types.STRING);
            put("xml", Types.OBJECT);
        }};
        descriptor.aggregations.add(new AggrFunctionInfo(Stats.STDEV, "stddev(#)", Types.dataFrameNumericTypes));
        descriptor.jdbcPropertiesTemplate = new ArrayList<Property>() {{
            // Application identification
            add(new Property(Property.STRING_TYPE, "applicationName",
                    "Name of the application shown in pg_stat_activity", new Prop()));

            // Timeouts
            add(new Property(Property.INT_TYPE, "connectTimeout",
                    "Connection timeout in seconds (0 = infinite, recommended > 0)", new Prop()));
            add(new Property(Property.INT_TYPE, "socketTimeout",
                    "Read timeout in seconds (0 = no timeout)", new Prop()));
            add(new Property(Property.INT_TYPE, "loginTimeout",
                    "Maximum time to wait for login in seconds", new Prop()));
            // Performance / batching
            add(new Property(Property.BOOL_TYPE, "reWriteBatchedInserts",
                    "Use optimized batched insert rewriting", new Prop()));
            add(new Property(Property.INT_TYPE, "prepareThreshold",
                    "Number of executes before using server-side prepared statements", new Prop()));
            // Misc
            add(new Property(Property.BOOL_TYPE, "tcpKeepAlive",
                    "Enable TCP keepalive", new Prop()));
        }};

    }

    @Override
    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = super.getProperties(conn);
        if (!conn.hasCustomConnectionString() && conn.ssl()) {
            properties.setProperty("ssl", "true");
            properties.setProperty("sslfactory", "org.postgresql.ssl.NonValidatingFactory");
        }
        if (!properties.containsKey("socketTimeout"))
            properties.setProperty("socketTimeout", "180");
        if (!properties.containsKey("tcpKeepAlive"))
            properties.setProperty("tcpKeepAlive", "true");
        return properties;
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:postgresql://" + conn.getServer() + port + "/" + conn.getDb();
    }

    @Override
    public String getSchemasSql(String db) {
        return "SELECT DISTINCT table_schema FROM information_schema.columns";
    }

    @Override
    public String getSchemaSql(String db, String schema, String table, boolean includeKeyInfo)
    {
        List<String> filters = new ArrayList<>();
        if (schema != null)
            filters.add("c.table_schema = '" + schema + "'");

        if (db != null && db.length() != 0)
            filters.add("c.table_catalog = '" + db + "'");

        if (table != null)
            filters.add("c.table_name = '" + table + "'");

        String whereClause = filters.isEmpty() ? "" : "WHERE " + String.join(" AND \n", filters);

        StringBuilder sql = new StringBuilder();

        sql.append("SELECT ")
                .append("c.table_catalog as table_catalog, ")
                .append("c.table_schema AS table_schema, ")
                .append("c.table_name AS table_name, ")
                .append("c.column_name AS column_name, ")
                .append("c.data_type AS data_type, ")
                .append("CASE t.table_type WHEN 'VIEW' THEN 1 ELSE 0 END AS is_view");


        if (includeKeyInfo) {
            sql.append(", ")
                    .append("CASE WHEN pk.column_name IS NOT NULL THEN 1 ELSE 0 END AS is_primary_key, ")
                    .append("CASE WHEN pk.column_name IS NOT NULL OR uq.column_name IS NOT NULL ")
                    .append("     THEN 1 ELSE 0 END AS is_unique");
        }

        sql.append(" FROM information_schema.columns c ")
                .append("JOIN information_schema.tables t ON t.table_name = c.table_name ")
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
                    .append("    AND pk.table_name = c.table_name ")
                    .append("    AND pk.column_name = c.column_name ");

            sql.append("LEFT JOIN ( ")
                    .append("    SELECT kcu.table_schema, kcu.table_name, kcu.column_name ")
                    .append("    FROM information_schema.table_constraints tc ")
                    .append("    JOIN information_schema.key_column_usage kcu ")
                    .append("      ON tc.constraint_name = kcu.constraint_name ")
                    .append("     AND tc.table_schema = kcu.table_schema ")
                    .append("    WHERE tc.constraint_type = 'UNIQUE' ")
                    .append(") uq ON uq.table_schema = c.table_schema ")
                    .append("    AND uq.table_name = c.table_name ")
                    .append("    AND uq.column_name = c.column_name ");
        }

        sql.append(" ")
                .append(whereClause)
                .append(" ORDER BY c.table_name, c.ordinal_position");

        return sql.toString();
    }

    @Override
    protected String getRegexQuery(String columnName, String regexExpression) {
        return String.format("%s ~ '%s'", columnName, regexExpression);
    }

    @Override
    protected void setUuid(PreparedStatement statement, int n, String value) throws SQLException {
        PGobject uuid = new PGobject();
        uuid.setType("uuid");
        uuid.setValue(value);
        statement.setObject(n, uuid);
    }

    @Override
    protected String getInQuery(PatternMatcher matcher, String names) {
        boolean isUuid = matcher.values != null && matcher.values.stream().allMatch((s) -> s instanceof String && UUID_REGEX.matcher((String) s).matches());
        return String.format("(%s%s %s (%s))", matcher.colName, isUuid ? "::uuid" : "", matcher.op, names);
    }

    @SuppressWarnings("unchecked")
    @Override
    protected int setArrayParamValue(PreparedStatement statement, int n, FuncParam param) throws SQLException {
        ArrayList<String> value = param.value == null ? null : (ArrayList<String>) param.value;
        if (value == null)
            statement.setNull(n, java.sql.Types.ARRAY);
        else {
            String type = value.isEmpty() || !value.stream().allMatch(s -> UUID_REGEX.matcher(s).matches())
                    ? "TEXT" : "UUID";
            Array array = statement.getConnection().createArrayOf(type, value.toArray());
            statement.setArray(n, array);
        }
        return 0;
    }

    @Override
    public String getCommentsQuery(DataConnection connection) {
        return "--input: string schema\n" +
                "SELECT \n" +
                "    n.nspname AS table_schema,\n" +
                "    'schema' AS object_type,\n" +
                "    NULL AS table_name,\n" +
                "    NULL AS column_name,\n" +
                "    obj_description(n.oid, 'pg_namespace') AS comment\n" +
                "FROM pg_namespace n\n" +
                "WHERE n.nspname = @schema\n" +
                "\n" +
                "UNION ALL\n" +
                "\n" +
                "SELECT\n" +
                "    n.nspname AS table_schema,\n" +
                "    CASE c.relkind\n" +
                "        WHEN 'r' THEN 'table'\n" +
                "        WHEN 'v' THEN 'view'\n" +
                "        WHEN 'm' THEN 'materialized_view'\n" +
                "        ELSE c.relkind::text\n" +
                "    END AS object_type,\n" +
                "    c.relname AS table_name,\n" +
                "    NULL AS column_name,\n" +
                "    obj_description(c.oid, 'pg_class') AS comment\n" +
                "FROM pg_class c\n" +
                "JOIN pg_namespace n ON n.oid = c.relnamespace\n" +
                "WHERE n.nspname = @schema\n" +
                "  AND c.relkind IN ('r','v','m')   -- r=table, v=view, m=mat view\n" +
                "\n" +
                "UNION ALL\n" +
                "\n" +
                "SELECT\n" +
                "    n.nspname AS table_schema,\n" +
                "    'column' AS object_type,\n" +
                "    c.relname AS table_name,\n" +
                "    a.attname AS column_name,\n" +
                "    col_description(c.oid, a.attnum) AS comment\n" +
                "FROM pg_class c\n" +
                "JOIN pg_attribute a ON a.attrelid = c.oid AND a.attnum > 0\n" +
                "JOIN pg_namespace n ON n.oid = c.relnamespace\n" +
                "WHERE n.nspname = @schema\n" +
                "  AND c.relkind IN ('r', 'v', 'm')\n" +
                "ORDER BY object_type, table_name, column_name;\n";
    }
}

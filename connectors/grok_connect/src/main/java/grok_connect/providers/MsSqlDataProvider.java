package grok_connect.providers;

import java.sql.*;
import java.util.*;

import grok_connect.managers.ColumnManager;
import grok_connect.managers.bool_column.MySqlMssqlBoolColumnManager;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataQuery;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncCall;
import grok_connect.connectors_info.FuncParam;
import grok_connect.resultset.DefaultResultSetManager;
import grok_connect.resultset.ResultSetManager;
import grok_connect.table_query.AggrFunctionInfo;
import grok_connect.table_query.Stats;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.GrokConnectUtil;
import grok_connect.utils.QueryCancelledByUser;
import serialization.DataFrame;
import serialization.Types;

public class MsSqlDataProvider extends JdbcDataProvider {

    public MsSqlDataProvider() {
        driverClassName = "com.microsoft.sqlserver.jdbc.SQLServerDriver";
        descriptor = new DataSource();
        descriptor.type = "MS SQL";
        descriptor.category = "Database";
        descriptor.description = "Query MS SQL database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(DbCredentials.getSsl());
        descriptor.credentialsTemplate = DbCredentials.getDbCredentialsTemplate();
        descriptor.canBrowseSchema = true;
        descriptor.supportCatalogs = true;
        descriptor.defaultSchema = "dbo";
        descriptor.limitAtEnd = false;

        descriptor.typesMap = new HashMap<String, String>() {{
            put("bigint", Types.BIG_INT);
            put("int", Types.INT);
            put("smallint", Types.INT);
            put("tinyint", Types.INT);
            put("#numeric.*", Types.FLOAT);
            put("#decimal.*", Types.FLOAT);
            put("money", Types.FLOAT);
            put("smallmoney", Types.FLOAT);
            put("float", Types.FLOAT);
            put("real", Types.FLOAT);
            put("bit", Types.BOOL);
            put("#.*char.*", Types.STRING);
            put("#.*varchar.*", Types.STRING);
            put("date", Types.DATE_TIME);
            put("datetimeoffset", Types.DATE_TIME);
            put("datetime2", Types.DATE_TIME);
            put("smalldatetime", Types.DATE_TIME);
            put("datetime", Types.DATE_TIME);
            put("time", Types.DATE_TIME);
            put("image", Types.OBJECT);
            put("#.*binary", Types.BLOB);
            put("#.*text.*", Types.OBJECT);
            put("geometry", Types.OBJECT);
            put("geography", Types.OBJECT);
            put("xml", Types.OBJECT);
        }};
        descriptor.aggregations.add(new AggrFunctionInfo(Stats.STDEV, "stdev(#)", Types.dataFrameNumericTypes));
    }

    @Override
    public String getConnectionString(DataConnection conn) {
        String connString = super.getConnectionString(conn);
        connString = connString.endsWith(";") ? connString : connString + ";";
        return connString;
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        String db = conn.getDb();
        String schema = conn.get("schema");
        db = (db == null || db.length() == 0) && schema != null ? schema : db;
        return "jdbc:sqlserver://" + conn.getServer() + port + ";databaseName=" + db +
                (conn.ssl() ? "integratedSecurity=true;encrypt=true;trustServerCertificate=true;" : "");
    }


    @Override
    public DataFrame getSchema(DataConnection connection, String schema, String table, boolean includeKeyInfo) throws QueryCancelledByUser,
            GrokConnectException {
        FuncCall queryRun = new FuncCall();
        queryRun.func = new DataQuery();
        queryRun.func.query = getSchemaSql(connection.getDb(), schema, table, includeKeyInfo);
        queryRun.func.connection = connection;

        return execute(queryRun);
    }

    @Override
    public String getSchemasSql(String db) {
        return "SELECT name as SCHEMA_NAME FROM sys.schemas ORDER BY name";
    }

    @Override
    public String getSchemaSql(String db, String schema, String table, boolean includeKeyInfo) {
        List<String> filters = new ArrayList<>();
        filters.add("o.type IN ('U', 'V', 'S')");

        if (schema != null)
            filters.add("s.name = '" + schema + "'");

        if (table != null)
            filters.add("o.name = '" + table + "'");

        String whereClause = "WHERE " + String.join(" AND \n", filters);

        StringBuilder sql = new StringBuilder();

        sql.append("SELECT ")
                .append("DB_NAME() AS table_catalog, ")
                .append("s.name AS table_schema, ")
                .append("o.name AS table_name, ")
                .append("c.name AS column_name, ")
                .append("ty.name AS data_type, ")
                .append("CASE WHEN o.type = 'V' THEN 1 ELSE 0 END AS is_view");

        if (includeKeyInfo) {
            sql.append(", ")
                    .append("CASE WHEN pk.column_id IS NOT NULL THEN 1 ELSE 0 END AS is_primary_key, ")
                    .append("CASE WHEN pk.column_id IS NOT NULL OR uq.column_id IS NOT NULL THEN 1 ELSE 0 END AS is_unique");
        }

        sql.append(" FROM sys.all_columns c ")
                .append("JOIN sys.all_objects o ON c.object_id = o.object_id ")
                .append("JOIN sys.schemas s ON o.schema_id = s.schema_id ")
                .append("JOIN sys.types ty ON c.user_type_id = ty.user_type_id ");

        if (includeKeyInfo) {
            sql.append("LEFT JOIN ( ")
                    .append("   SELECT ic.object_id, ic.column_id ")
                    .append("   FROM sys.index_columns ic ")
                    .append("   JOIN sys.indexes i ON ic.object_id = i.object_id AND ic.index_id = i.index_id ")
                    .append("   WHERE i.is_primary_key = 1 ")
                    .append(") pk ON pk.object_id = c.object_id AND pk.column_id = c.column_id ");

            sql.append("LEFT JOIN ( ")
                    .append("   SELECT ic.object_id, ic.column_id ")
                    .append("   FROM sys.index_columns ic ")
                    .append("   JOIN sys.indexes i ON ic.object_id = i.object_id AND ic.index_id = i.index_id ")
                    .append("   WHERE i.is_unique = 1 ")
                    .append(") uq ON uq.object_id = c.object_id AND uq.column_id = c.column_id ");
        }

        sql.append(" ").append(whereClause);
        sql.append(" ORDER BY o.name, c.column_id");

        return sql.toString();
    }

    @Override
    protected void appendQueryParam(DataQuery dataQuery, String paramName, StringBuilder queryBuffer) {
        FuncParam param = dataQuery.getParam(paramName);
        if (param.propertyType.equals("list")) {
            queryBuffer.append("SELECT value FROM STRING_SPLIT(?, ',')");
        } else {
            queryBuffer.append("?");
        }
    }

    @Override
    protected int setArrayParamValue(PreparedStatement statement, int n, FuncParam param) throws SQLException {
        @SuppressWarnings("unchecked")
        List<String> list = ((ArrayList<String>) param.value);
        String values = String.join(",", list);
        statement.setObject(n, values);
        return 0;
    }

    @Override
    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = defaultConnectionProperties(conn);
        properties.setProperty("socketTimeout", "180000");
        return properties;
    }

    @Override
    public String limitToSql(String query, Integer limit) {
        return String.format("%stop %s", query, limit);
    }

    @Override
    public ResultSetManager getResultSetManager() {
        Map<String, ColumnManager<?>> defaultManagersMap = DefaultResultSetManager.getDefaultManagersMap();
        defaultManagersMap.put(Types.BOOL, new MySqlMssqlBoolColumnManager());
        return DefaultResultSetManager.fromManagersMap(defaultManagersMap);
    }


    @Override
    public DataFrame getForeignKeys(DataConnection conn, String schema) throws GrokConnectException, QueryCancelledByUser {
        FuncCall queryRun = new FuncCall();
        queryRun.func = new DataQuery();
        queryRun.func.query = "SELECT sch.name AS [table_schema], \n" +
                "    obj.name AS [constraint_name],\n" +
                "        tab1.name AS [table_name],\n" +
                "        col1.name AS [column_name],\n" +
                "        tab2.name AS [foreign_table_name],\n" +
                "        col2.name AS [foreign_column_name]\n" +
                "    FROM sys.foreign_key_columns fkc\n" +
                "    INNER JOIN sys.objects obj\n" +
                "        ON obj.object_id = fkc.constraint_object_id\n" +
                "    INNER JOIN sys.tables tab1\n" +
                "        ON tab1.object_id = fkc.parent_object_id\n" +
                "    INNER JOIN sys.schemas sch\n" +
                "        ON tab1.schema_id = sch.schema_id\n" +
                "    INNER JOIN sys.columns col1\n" +
                "        ON col1.column_id = parent_column_id AND col1.object_id = tab1.object_id\n" +
                "    INNER JOIN sys.tables tab2\n" +
                "        ON tab2.object_id = fkc.referenced_object_id\n" +
                "    INNER JOIN sys.columns col2\n" +
                "        ON col2.column_id = referenced_column_id AND col2.object_id = tab2.object_id";
        if (GrokConnectUtil.isNotEmpty(schema))
            queryRun.func.query += String.format("\n\t\tWHERE sch.name = '%s'", schema);
        queryRun.func.connection = conn;
        return execute(queryRun);
    }

    @Override
    public String getCommentsQuery(DataConnection connection) throws GrokConnectException {
        return "--input: string schema\n" +
                "SELECT\n" +
                "    s.name AS table_schema,\n" +
                "    'schema' AS object_type,\n" +
                "    NULL AS table_name,\n" +
                "    NULL AS column_name,\n" +
                "    CAST(ep.value AS NVARCHAR(MAX)) AS comment\n" +
                "FROM sys.schemas s\n" +
                "OUTER APPLY fn_listextendedproperty(\n" +
                "        'MS_Description',\n" +
                "        'schema', s.name,\n" +
                "        NULL, NULL,\n" +
                "        NULL, NULL\n" +
                "     ) ep\n" +
                "WHERE s.name = @schema\n" +
                "\n" +
                "UNION ALL\n" +
                "\n" +
                "SELECT\n" +
                "    s.name AS table_schema,\n" +
                "    'table' AS object_type,\n" +
                "    t.name AS table_name,\n" +
                "    NULL AS column_name,\n" +
                "    CAST(ep.value AS NVARCHAR(MAX)) AS comment\n" +
                "FROM sys.tables t\n" +
                "JOIN sys.schemas s ON s.schema_id = t.schema_id\n" +
                "OUTER APPLY fn_listextendedproperty(\n" +
                "        'MS_Description',\n" +
                "        'schema', s.name,\n" +
                "        'table', t.name,\n" +
                "        NULL, NULL\n" +
                "     ) ep\n" +
                "WHERE s.name = @schema\n" +
                "\n" +
                "UNION ALL\n" +
                "\n" +
                "SELECT\n" +
                "    s.name AS table_schema,\n" +
                "    'view' AS object_type,\n" +
                "    v.name AS table_name,\n" +
                "    NULL AS column_name,\n" +
                "    CAST(ep.value AS NVARCHAR(MAX)) AS comment\n" +
                "FROM sys.views v\n" +
                "JOIN sys.schemas s ON s.schema_id = v.schema_id\n" +
                "OUTER APPLY fn_listextendedproperty(\n" +
                "        'MS_Description',\n" +
                "        'schema', s.name,\n" +
                "        'view', v.name,\n" +
                "        NULL, NULL\n" +
                "     ) ep\n" +
                "WHERE s.name = @schema\n" +
                "\n" +
                "UNION ALL\n" +
                "\n" +
                "SELECT\n" +
                "    s.name AS table_schema,\n" +
                "    'column' AS object_type,\n" +
                "    obj.name AS table_name,\n" +
                "    c.name AS column_name,\n" +
                "    CAST(ep.value AS NVARCHAR(MAX)) AS comment\n" +
                "FROM sys.columns c\n" +
                "JOIN sys.objects obj ON obj.object_id = c.object_id\n" +
                "JOIN sys.schemas s ON s.schema_id = obj.schema_id\n" +
                "OUTER APPLY fn_listextendedproperty(\n" +
                "        'MS_Description',\n" +
                "        'schema', s.name,\n" +
                "        CASE WHEN obj.type = 'V' THEN 'view' ELSE 'table' END,\n" +
                "        obj.name,\n" +
                "        'column', c.name\n" +
                "     ) ep\n" +
                "WHERE s.name = @schema\n" +
                "  AND obj.type IN ('U', 'V')   -- U = tables, V = views\n" +
                "\n" +
                "ORDER BY object_type, table_name, column_name;\n";
    }
}

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
import grok_connect.utils.Property;
import grok_connect.utils.QueryCancelledByUser;
import serialization.Types;
import serialization.DataFrame;

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
    public DataFrame getSchema(DataConnection connection, String schema, String table) throws QueryCancelledByUser,
            GrokConnectException {
        FuncCall queryRun = new FuncCall();
        queryRun.func = new DataQuery();
        String db = connection.getDb();
        queryRun.func.query = (db != null && db.length() != 0)
                ? getSchemaSql(db, schema, table) : getSchemaSql(schema, null, table);
        queryRun.func.connection = connection;

        return execute(queryRun);
    }

    @Override
    public String getSchemasSql(String db) {
        return "SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA";
    }

    @Override
    public String getSchemaSql(String db, String schema, String table) {
        List<String> filters = new ArrayList<String>() {{
            add("c.table_schema = '" + ((schema != null) ? schema : descriptor.defaultSchema) + "'");
        }};

        if (db != null && db.length() != 0)
            filters.add("c.table_catalog = '" + db + "'");

        if (table!= null)
            filters.add("c.table_name = '" + table + "'");

        String whereClause = "WHERE " + String.join(" AND \n", filters);

        return "SELECT c.table_schema as table_schema, c.table_name as table_name, c.column_name as column_name, "
                + "c.data_type as data_type, "
                + "case t.table_type when 'VIEW' then 1 else 0 end as is_view FROM information_schema.columns c "
                + "JOIN information_schema.tables t ON t.table_name = c.table_name AND t.table_schema = c.table_schema AND t.table_catalog = c.table_catalog " + whereClause +
                " ORDER BY c.table_name";
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
}

package grok_connect.providers;

import java.io.*;
import java.sql.*;
import java.text.*;
import java.util.*;
import serialization.Types;
import serialization.DataFrame;
import grok_connect.utils.*;
import grok_connect.table_query.*;
import grok_connect.connectors_info.*;


public class MsSqlDataProvider extends JdbcDataProvider {
    public MsSqlDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "MS SQL";
        descriptor.category = "Database";
        descriptor.description = "Query MS SQL database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.canBrowseSchema = true;
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
            put("#.char.*", Types.STRING);
            put("#.varchar.*", Types.STRING);
            put("#.text.*", Types.STRING);
            put("date", Types.DATE_TIME);
            put("datetimeoffset", Types.DATE_TIME);
            put("datetime2", Types.DATE_TIME);
            put("smalldatetime", Types.DATE_TIME);
            put("datetime", Types.DATE_TIME);
            put("time", Types.DATE_TIME);
        }};
        descriptor.aggregations = new ArrayList<AggrFunctionInfo>() {{
            add(new AggrFunctionInfo(Stats.AVG, "avg(#)", Types.dataFrameNumericTypes));
            add(new AggrFunctionInfo(Stats.MIN, "min(#)", Types.dataFrameNumericTypes));
            add(new AggrFunctionInfo(Stats.MAX, "max(#)", Types.dataFrameNumericTypes));
            add(new AggrFunctionInfo(Stats.SUM, "sum(#)", Types.dataFrameNumericTypes));
            add(new AggrFunctionInfo(Stats.STDEV, "stdev(#)", Types.dataFrameNumericTypes));
            add(new AggrFunctionInfo(Stats.TOTAL_COUNT, "count(*)", Types.dataFrameColumnTypes));
            add(new AggrFunctionInfo(Stats.VALUE_COUNT, "count(#)", Types.dataFrameColumnTypes));
            add(new AggrFunctionInfo(Stats.MISSING_VALUE_COUNT, "count(*) - count(#)", Types.dataFrameColumnTypes));
        }};
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("com.microsoft.sqlserver.jdbc.SQLServerDriver");
        String connString = getConnectionString(conn);
        connString = connString.endsWith(";") ? connString : connString + ";";
        connString += "user=" + conn.credentials.getLogin() + ";password=" + conn.credentials.getPassword();
        return DriverManager.getConnection(connString);
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        String db = conn.getDb();
        String schema = conn.get("schema");
        db = (db == null || db.length() == 0) && schema != null ? schema : db;
        return "jdbc:sqlserver://" + conn.getServer() + port + ";databaseName=" + db +
                (conn.ssl() ? "integratedSecurity=true;encrypt=true;trustServerCertificate=true;" : "");
    }

    public DataFrame getSchema(DataConnection connection, String schema, String table)
            throws ClassNotFoundException, SQLException, ParseException, IOException {
        FuncCall queryRun = new FuncCall();
        queryRun.func = new DataQuery();
        queryRun.func.query = getSchemaSql(schema, null, table);
        queryRun.func.connection = connection;

        return execute(queryRun);
    }

    public String getSchemasSql(String db) {
        return "SELECT DISTINCT name FROM sys.databases";
    }

    public String getSchemaSql(String db, String schema, String table)
    {
        List<String> filters = new ArrayList<String>() {{
            add("TABLE_SCHEMA = '" + ((schema != null) ? schema : "dbo") + "'");
        }};

        if (db != null && db.length() != 0)
            filters.add("TABLE_CATALOG = '" + db + "'");

        if (table!= null)
            filters.add("TABLE_NAME = '" + table + "'");

        String whereClause = "WHERE " + String.join(" AND \n", filters);

        return "SELECT TABLE_SCHEMA, TABLE_NAME, COLUMN_NAME, DATA_TYPE " +
                "FROM INFORMATION_SCHEMA.COLUMNS " + whereClause + " ORDER BY TABLE_NAME";
    }

    public String limitToSql(String query, Integer limit) {
        return query + "top " + limit.toString() + " ";
    }
}

package grok_connect.providers;

import java.sql.*;
import java.util.*;
import serialization.Types;
import grok_connect.utils.*;
import grok_connect.table_query.*;
import grok_connect.connectors_info.*;


public class MsSqlDataProvider extends JdbcDataProvider {
    public MsSqlDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "MS SQL";
        descriptor.category = "Database";
        descriptor.description = "Query MS SQL database";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
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
        return DriverManager.getConnection(getConnectionString(conn));
    }

    public String getConnectionString(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        boolean ssl = conn.parameters.containsKey(DbCredentials.SSL) && (boolean)conn.parameters.get(DbCredentials.SSL);
        return "jdbc:sqlserver://" + conn.getServer() + port + ";databaseName=" + conn.getDb() +
                ";user=" + conn.credentials.getLogin() + ";password=" + conn.credentials.getPassword() + ";" +
                (ssl ? "integratedSecurity=true;encrypt=true;trustServerCertificate=true;" : "");
    }

    public String getSchemasSql(String db) {
        return "SELECT DISTINCT TABLE_SCHEMA FROM INFORMATION_SCHEMA.COLUMNS ORDER BY TABLE_SCHEMA";
    }

    public String getSchemaSql(String db, String schema, String table)
    {
        List<String> filters = new ArrayList<>();
        filters.add("(TABLE_SCHEMA = '" + ((schema != null) ? schema : "dbo") + "')");
        if (db != null)
            filters.add("(TABLE_CATALOG = '" + db + "')");
        if (table!= null)
            filters.add("(TABLE_NAME = '" + table + "')");

        String whereClause = String.join(" AND \n", filters);

        return "SELECT TABLE_SCHEMA, TABLE_NAME, COLUMN_NAME, DATA_TYPE " +
                "FROM INFORMATION_SCHEMA.COLUMNS WHERE " + whereClause + " ORDER BY TABLE_NAME";
    }

    public String limitToSql(String query, Integer limit) {
        return query + "top " + limit.toString() + " ";
    }
}

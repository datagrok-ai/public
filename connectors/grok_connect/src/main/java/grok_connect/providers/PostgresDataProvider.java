package grok_connect.providers;

import java.sql.*;
import java.util.*;
import serialization.Types;
import grok_connect.utils.*;
import grok_connect.table_query.*;
import grok_connect.connectors_info.*;


public class PostgresDataProvider extends JdbcDataProvider {
    private static final String SYS_SCHEMAS_FILTER =
            "table_schema != 'pg_catalog' AND table_schema != 'information_schema'";

    public PostgresDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "PostgresNet";
        descriptor.description = "Query PostgresNet database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.nameBrackets = "\"";

        descriptor.canBrowseSchema = true;
        descriptor.typesMap = new HashMap<String, String>() {{
            put("smallint", Types.INT);
            put("int", Types.INT);
            put("bigint", Types.BIG_INT);
            put("real", Types.FLOAT);
            put("double precision", Types.FLOAT);
            put("numeric", Types.FLOAT);
            put("#character.*", Types.STRING);
            put("#varchar.*", Types.STRING);
            put("text", Types.STRING);
        }};
        descriptor.aggregations = new ArrayList<AggrFunctionInfo>() {{
            add(new AggrFunctionInfo(Stats.AVG, "avg(#)", Types.dataFrameNumericTypes));
            add(new AggrFunctionInfo(Stats.MIN, "min(#)", Types.dataFrameNumericTypes));
            add(new AggrFunctionInfo(Stats.MAX, "max(#)", Types.dataFrameNumericTypes));
            add(new AggrFunctionInfo(Stats.SUM, "sum(#)", Types.dataFrameNumericTypes));
            add(new AggrFunctionInfo(Stats.STDEV, "stddev(#)", Types.dataFrameNumericTypes));
            add(new AggrFunctionInfo(Stats.TOTAL_COUNT, "count(*)", Types.dataFrameColumnTypes));
            add(new AggrFunctionInfo(Stats.VALUE_COUNT, "count(#)", Types.dataFrameColumnTypes));
            add(new AggrFunctionInfo(Stats.MISSING_VALUE_COUNT, "count(*) - count(#)", Types.dataFrameColumnTypes));
        }};
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("org.postgresql.Driver");
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString() && conn.ssl())
            properties.setProperty("ssl", "true");
        return DriverManager.getConnection(getConnectionString(conn), properties);
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:postgresql://" + conn.getServer() + port + "/" + conn.getDb();
    }

    public String getSchemasSql(String db) {
        return "SELECT DISTINCT table_schema FROM information_schema.columns " +
                "WHERE table_catalog = '" + db + "' AND " + SYS_SCHEMAS_FILTER + " ORDER BY table_schema";
    }

    public String getSchemaSql(String db, String schema, String table)
    {
        List<String> filters = new ArrayList<String>() {{
            add("table_schema = '" + ((schema != null) ? schema : "public") + "'");
            add("table_catalog = '" + db + "'");
        }};

        if (table != null)
            filters.add("table_name = '" + table + "'");

        String whereClause = String.join(" and \n", filters);

        return "SELECT table_schema, table_name, column_name, data_type " +
                "FROM information_schema.columns WHERE " + SYS_SCHEMAS_FILTER + " AND " + whereClause +
                " ORDER BY table_name";
    }
}

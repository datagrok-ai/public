package grok_connect.providers;

import java.sql.*;
import java.util.*;

import serialization.Types;
import grok_connect.table_query.*;
import grok_connect.connectors_info.*;


public class OracleDataProvider extends JdbcDataProvider {
    public OracleDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Oracle";
        descriptor.description = "Query Oracle database";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
        descriptor.canBrowseSchema = true;

        descriptor.typesMap = new HashMap<String, String>() {{
            put("long", Types.INT);
            put("#float.*", Types.FLOAT);
            put("#number.*", Types.FLOAT);
            put("binary_float", Types.FLOAT);
            put("binary_double", Types.FLOAT);
            put("#.char.*", Types.STRING);
            put("#.varchar.*", Types.STRING);
            put("#.clob.*", Types.STRING);
            put("date", Types.DATE_TIME);
            put("timestamp", Types.DATE_TIME);
        }};
        descriptor.aggregations = new ArrayList<AggrFunctionInfo>() {{
            add(new AggrFunctionInfo(Stats.AVG, "avg(#)", Types.dataFrameNumericTypes));
            add(new AggrFunctionInfo(Stats.MIN, "min(#)", Types.dataFrameNumericTypes));
            add(new AggrFunctionInfo(Stats.MAX, "max(#)", Types.dataFrameNumericTypes));
            add(new AggrFunctionInfo(Stats.SUM, "sum(#)", Types.dataFrameNumericTypes));
            add(new AggrFunctionInfo(Stats.STDEV, "stddev(#)", Types.dataFrameNumericTypes));
            add(new AggrFunctionInfo(Stats.MED, "median(#)", Types.dataFrameNumericTypes));
            add(new AggrFunctionInfo(Stats.TOTAL_COUNT, "count(*)", Types.dataFrameColumnTypes));
            add(new AggrFunctionInfo(Stats.VALUE_COUNT, "count(#)", Types.dataFrameColumnTypes));
            add(new AggrFunctionInfo(Stats.MISSING_VALUE_COUNT, "count(*) - count(#)", Types.dataFrameColumnTypes));
        }};
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("oracle.jdbc.driver.OracleDriver");
        return DriverManager.getConnection(getConnectionString(conn), conn.getLogin(), conn.getPassword());
    }

    public String getConnectionString(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:oracle:thin:@//" + conn.getServer() + port + "/" + conn.getDb();
    }

    public String getSchemaSql(String db, String schema, String table)
    {
        List<String> filters = new ArrayList<>();
        if (table!= null)
            filters.add("(TABLE_NAME = '" + table + "')");

        String whereClause = " WHERE " + String.join(" AND \n", filters);

        return "SELECT TABLE_NAME, COLUMN_NAME, DATA_TYPE " +
                "FROM USER_TAB_COLUMNS" + ((filters.size() > 0) ? whereClause : "");
    }

    public String limitToSql(String query, Integer limit) {
        return "select * from (\n" + query + "\n) where ROWNUM <= " + limit.toString();
    }
}

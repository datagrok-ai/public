package grok_connect.providers;

import java.sql.*;
import java.util.*;

import grok_connect.utils.*;
import serialization.Types;
import grok_connect.table_query.*;
import grok_connect.connectors_info.*;


public class OracleDataProvider extends JdbcDataProvider {
    private static final String SYS_SCHEMAS_FILTER =
            "OWNER != 'SYSTEM' AND OWNER != 'CTXSYS' AND OWNER != 'MDSYS' " +
            "AND OWNER != 'XDB' AND OWNER != 'APEX_040000' AND OWNER != 'SYS' " +
            "AND OWNER != 'WMSYS' AND OWNER != 'EXFSYS' AND OWNER != 'ORDSYS' " +
            "AND OWNER != 'ORDDATA'";

    public OracleDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Oracle";
        descriptor.description = "Query Oracle database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.canBrowseSchema = true;
        descriptor.nameBrackets = "\"";

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
        descriptor.aggregations.add(new AggrFunctionInfo(Stats.STDEV, "stddev(#)", Types.dataFrameNumericTypes));
        descriptor.aggregations.add(new AggrFunctionInfo(Stats.STDEV, "median(#)", Types.dataFrameNumericTypes));
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("oracle.jdbc.driver.OracleDriver");
        System.getProperties().setProperty("oracle.jdbc.J2EE13Compliant", "true");
        return DriverManager.getConnection(getConnectionString(conn), conn.credentials.getLogin(), conn.credentials.getPassword());
    }

    public String getConnectionStringImpl(DataConnection conn) {
        conn.getPort();
        return "jdbc:oracle:thin:@(DESCRIPTION=" +
                "(ADDRESS=" +
                    "(PROTOCOL=" + (conn.ssl() ? "tcps" : "tcp") + ")" +
                    "(HOST=" + conn.getServer() + ")" +
                    "(PORT=" + conn.getPort() + "))" +
                "(CONNECT_DATA=(SERVICE_NAME=" + conn.getDb() + ")))";
    }

    public String getSchemasSql(String db) {
        return "SELECT OWNER as TABLE_SCHEMA FROM ALL_TAB_COLUMNS WHERE " + SYS_SCHEMAS_FILTER +
                " GROUP BY OWNER ORDER BY OWNER";
    }

    public String getSchemaSql(String db, String schema, String table) {
        String whereClause = "WHERE " + SYS_SCHEMAS_FILTER;

        if (table != null)
            whereClause = whereClause + " AND (TABLE_NAME = '" + table + "')";
        if (schema != null)
            whereClause = whereClause + " AND (OWNER = '" + schema + "')";

        return "SELECT OWNER as TABLE_SCHEMA, TABLE_NAME, COLUMN_NAME, DATA_TYPE FROM ALL_TAB_COLUMNS " + whereClause +
                " ORDER BY TABLE_NAME";
    }

    public String limitToSql(String query, Integer limit) {
        return "select * from (\n" + query + "\n) where ROWNUM <= " + limit.toString();
    }

    public String addBrackets(String name) {
        String brackets = descriptor.nameBrackets;
        return name.startsWith(brackets.substring(0, 1)) ? name :
                brackets.substring(0, 1) + name + brackets.substring(brackets.length() - 1, brackets.length());
    }

    private boolean isOracleFloatNumber(String typeName, int precision, int scale) {
        // https://markhoxey.wordpress.com/2016/05/31/maximum-number-precision/ ==>  Precision >= 38
        // https://stackoverflow.com/questions/29537292/why-can-number-type-in-oracle-scale-up-to-127 ==> scale >= 127
        return typeName.equalsIgnoreCase("number") && (scale == 0 || scale >= 127) && (precision <= 0);
    }

    @Override
    protected boolean isInteger(int type, String typeName, int precision, int scale) {
        // https://docs.oracle.com/cd/E11882_01/server.112/e41084/sql_elements001.htm#sthref119
        // The absence of precision and scale designators specifies the maximum range and precision for an Oracle number.
        // We shall ignore the case where type == java.sql.Types. ... value is identified incorrectly
        if (isOracleFloatNumber(typeName, precision, scale)) return false;
        return super.isInteger(type, typeName, precision, scale);
    }

    @Override
    protected boolean isFloat(int type, String typeName, int precision, int scale) {
        return super.isFloat(type, typeName, precision, scale) || isOracleFloatNumber(typeName, precision, scale);
    }
}

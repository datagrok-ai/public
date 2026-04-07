package grok_connect.providers;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.*;
import java.util.stream.Collectors;
import grok_connect.managers.ColumnManager;
import grok_connect.managers.bigint_column.OracleBigIntColumnManager;
import grok_connect.managers.integer_column.OracleSnowflakeIntColumnManager;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataQuery;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncParam;
import grok_connect.resultset.DefaultResultSetManager;
import grok_connect.resultset.OracleResultSetManager;
import grok_connect.resultset.ResultSetManager;
import grok_connect.table_query.AggrFunctionInfo;
import grok_connect.table_query.Stats;
import grok_connect.utils.Prop;
import grok_connect.utils.Property;
import oracle.jdbc.OracleResultSet;
import oracle.jdbc.driver.OracleConnection;
import oracle.sql.json.OracleJsonObject;
import serialization.Types;

public class OracleDataProvider extends JdbcDataProvider {
    private static final String SYS_SCHEMAS_FILTER =
            "COL.OWNER != 'SYSTEM' AND COL.OWNER != 'CTXSYS' AND COL.OWNER != 'MDSYS' " +
            "AND COL.OWNER != 'XDB' AND COL.OWNER != 'APEX_040000' AND COL.OWNER != 'SYS' " +
            "AND COL.OWNER != 'WMSYS' AND COL.OWNER != 'EXFSYS' AND COL.OWNER != 'ORDSYS' " +
            "AND COL.OWNER != 'ORDDATA'";

    public OracleDataProvider() {
        driverClassName = "oracle.jdbc.OracleDriver";
        descriptor = new DataSource();
        descriptor.type = "Oracle";
        descriptor.description = "Query Oracle database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(DbCredentials.getSsl());
        descriptor.credentialsTemplate = DbCredentials.getDbCredentialsTemplate();
        descriptor.canBrowseSchema = true;
        descriptor.nameBrackets = "\"";

        descriptor.typesMap = new HashMap<String, String>() {{
            put("long", Types.INT);
            put("float", Types.FLOAT);
            put("#number\\([1-9], 0\\)", Types.INT);
            put("#number\\((3[0-8]|[1-2][0-9]), 0\\)", Types.BIG_INT);
            put("#number\\((3[0-8]|[1-2][0-9]|[1-9]), -?[^0]+\\)", Types.FLOAT);
            put("binary_float", Types.FLOAT);
            put("binary_double", Types.FLOAT);
            put("#.*char.*", Types.STRING);
            put("#.*varchar.*", Types.STRING);
            put("date", Types.DATE_TIME);
            put("#timestamp.*", Types.DATE_TIME);
            put("#interval.*", Types.DATE_TIME);
            put("json", Types.OBJECT);
            put("#.*clob.*", Types.BLOB);
            put("blob", Types.BLOB);
            put("uritype", Types.OBJECT);
            put("mem_type", Types.OBJECT);
            put("xmltype", Types.OBJECT);
        }};
        descriptor.aggregations.add(new AggrFunctionInfo(Stats.STDEV, "stddev(#)", Types.dataFrameNumericTypes));
        descriptor.aggregations.add(new AggrFunctionInfo(Stats.STDEV, "median(#)", Types.dataFrameNumericTypes));
        descriptor.jdbcPropertiesTemplate = new ArrayList<Property>() {{
            add(new Property(Property.INT_TYPE, OracleConnection.CONNECTION_PROPERTY_LOGIN_TIMEOUT,
                    "Maximum time to wait for login in seconds", new Prop(), "loginTimeout"));
            add(new Property(Property.INT_TYPE, OracleConnection.CONNECTION_PROPERTY_THIN_NET_CONNECT_TIMEOUT,
                    "Connection timeout in milliseconds", new Prop(), "connectTimeout"));
            add(new Property(Property.INT_TYPE, OracleConnection.CONNECTION_PROPERTY_THIN_READ_TIMEOUT,
                    "Socket read timeout in milliseconds", new Prop(), "readTimeout"));

            add(new Property(Property.BOOL_TYPE, OracleConnection.CONNECTION_PROPERTY_REPORT_REMARKS,
                    "Enable retrieval of table/column remarks via DatabaseMetaData", new Prop(), "remarksReporting"));

            add(new Property(Property.INT_TYPE, OracleConnection.CONNECTION_PROPERTY_IMPLICIT_STATEMENT_CACHE_SIZE,
                    "Number of statements in the implicit cache", new Prop(), "implicitStatementCacheSize"));
            add(new Property(Property.INT_TYPE, OracleConnection.CONNECTION_PROPERTY_MAX_CACHED_BUFFER_SIZE,
                    "Maximum cached buffer size in KB", new Prop(), "maxCachedBufferSize"));
            add(new Property(Property.BOOL_TYPE, OracleConnection.CONNECTION_PROPERTY_USE_FETCH_SIZE_WITH_LONG_COLUMN,
                    "Apply fetch size to LONG and LONG RAW columns", new Prop()));

            add(new Property(Property.BOOL_TYPE, OracleConnection.CONNECTION_PROPERTY_DEFAULTNCHAR,
                    "Use NCHAR semantics for all character data", new Prop()));
        }};

    }

    @Override
    public void prepareProvider() {
        System.getProperties().setProperty("oracle.jdbc.J2EE13Compliant", "true");
    }

    @Override
    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = super.getProperties(conn);
        if (!properties.containsKey(OracleConnection.CONNECTION_PROPERTY_THIN_NET_CONNECT_TIMEOUT))
            properties.setProperty(OracleConnection.CONNECTION_PROPERTY_THIN_NET_CONNECT_TIMEOUT, "180000");
        return properties;
    }

    @Override
    protected String getRegexQuery(String columnName, String regexExpression) {
        return String.format("REGEXP_LIKE (%s, '%s', 'i')", columnName, regexExpression);
    }

    @Override
    protected Object getObjectFromResultSet(ResultSet resultSet, int c) {
        try {
            if (resultSet.getMetaData().getColumnTypeName(c).equals("JSON")) {
                return resultSet.unwrap(OracleResultSet.class).getObject(c, OracleJsonObject.class);
            }
            return resultSet.getObject(c);
        } catch (SQLException e) {
            throw new RuntimeException("Something went wrong when getting object from result set");
        }
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
    public String getConnectionStringImpl(DataConnection conn) {
        conn.getPort();
        return "jdbc:oracle:thin:@(DESCRIPTION=" +
                "(ADDRESS=" +
                    "(PROTOCOL=" + (conn.ssl() ? "tcps" : "tcp") + ")" +
                    "(HOST=" + conn.getServer() + ")" +
                    "(PORT=" + conn.getPort() + "))" +
                "(CONNECT_DATA=(SERVICE_NAME=" + conn.getDb() + ")))";
    }

    @Override
    public String getSchemasSql(String db) {
        return "SELECT COL.OWNER as TABLE_SCHEMA FROM ALL_TAB_COLUMNS COL WHERE " + SYS_SCHEMAS_FILTER +
                " GROUP BY COL.OWNER ORDER BY COL.OWNER";
    }

    @Override
    public String getSchemaSql(String db, String schema, String table, boolean includeKeyInfo) {
        String whereClause = "WHERE " + SYS_SCHEMAS_FILTER;

        if (table != null)
            whereClause = whereClause + " AND (COL.TABLE_NAME = '" + table + "')";
        if (schema != null)
            whereClause = whereClause + " AND (COL.OWNER = '" + schema + "')";

        StringBuilder sql = new StringBuilder();

        sql.append(
                "SELECT " +
                        "  COL.OWNER AS TABLE_SCHEMA, " +
                        "  COL.TABLE_NAME AS TABLE_NAME, " +
                        "  COL.COLUMN_NAME AS COLUMN_NAME, " +
                        "  CASE " +
                        "    WHEN DATA_PRECISION IS NOT NULL AND DATA_SCALE IS NOT NULL THEN " +
                        "      COL.DATA_TYPE || '(' || DATA_PRECISION || ', ' || DATA_SCALE || ')' " +
                        "    ELSE COL.DATA_TYPE " +
                        "  END AS DATA_TYPE, " +
                        "  CASE WHEN O.OBJECT_TYPE = 'VIEW' THEN 1 ELSE 0 END AS IS_VIEW"
        );

        if (includeKeyInfo) {
            sql.append(
                    ", CASE WHEN pk.column_name IS NOT NULL THEN 1 ELSE 0 END AS is_primary_key " +
                            ", CASE WHEN pk.column_name IS NOT NULL OR uq.column_name IS NOT NULL " +
                            "       THEN 1 ELSE 0 END AS is_unique"
            );
        }

        sql.append(
                " FROM ALL_TAB_COLUMNS COL " +
                        " JOIN ALL_OBJECTS O ON O.OBJECT_NAME = COL.TABLE_NAME AND O.OWNER = COL.OWNER "
        );

        if (includeKeyInfo) {
            sql.append(
                    " LEFT JOIN ( " +
                            "   SELECT acc.owner, acc.table_name, acc.column_name " +
                            "   FROM ALL_CONSTRAINTS ac " +
                            "   JOIN ALL_CONS_COLUMNS acc ON ac.constraint_name = acc.constraint_name AND ac.owner = acc.owner " +
                            "   WHERE ac.constraint_type = 'P' " + // Primary key
                            " ) pk ON pk.owner = COL.OWNER AND pk.table_name = COL.TABLE_NAME AND pk.column_name = COL.COLUMN_NAME "
            );

            sql.append(
                    " LEFT JOIN ( " +
                            "   SELECT acc.owner, acc.table_name, acc.column_name " +
                            "   FROM ALL_CONSTRAINTS ac " +
                            "   JOIN ALL_CONS_COLUMNS acc ON ac.constraint_name = acc.constraint_name AND ac.owner = acc.owner " +
                            "   WHERE ac.constraint_type = 'U' " + // Unique
                            " ) uq ON uq.owner = COL.OWNER AND uq.table_name = COL.TABLE_NAME AND uq.column_name = COL.COLUMN_NAME "
            );
        }

        sql.append(" ")
                .append(whereClause)
                .append(" ORDER BY TABLE_NAME");

        return sql.toString();
    }


    @Override
    public String limitToSql(String query, Integer limit) {
        return String.format("SELECT * FROM (%s%s%s) WHERE ROWNUM <= %s", System.lineSeparator(),
                query, System.lineSeparator(), limit);
    }

    @Override
    public ResultSetManager getResultSetManager() {
        Map<String, ColumnManager<?>> defaultManagersMap = DefaultResultSetManager.getDefaultManagersMap();
        defaultManagersMap.put(Types.INT, new OracleSnowflakeIntColumnManager());
        defaultManagersMap.put(Types.BIG_INT, new OracleBigIntColumnManager());
        return new OracleResultSetManager(defaultManagersMap.values());
    }

    @Override
    public String getCommentsQuery(DataConnection connection) {
        return "--input: string schema\n" +
                "SELECT\n" +
                "    tc.owner AS table_schema,\n" +
                "    tc.table_name AS object_name,\n" +
                "    CASE \n" +
                "        WHEN t.table_type = 'VIEW' THEN 'view'\n" +
                "        ELSE 'table'\n" +
                "    END AS object_type,\n" +
                "    NULL AS column_name,\n" +
                "    tc.comments AS comment\n" +
                "FROM all_tab_comments tc\n" +
                "JOIN all_objects t\n" +
                "  ON t.owner = tc.owner\n" +
                " AND t.object_name = tc.table_name\n" +
                "WHERE tc.owner = @schema\n" +
                "\n" +
                "UNION ALL\n" +
                "\n" +
                "-- COLUMN comments\n" +
                "SELECT\n" +
                "    cc.owner AS table_schema,\n" +
                "    cc.table_name AS object_name,\n" +
                "    'column' AS object_type,\n" +
                "    cc.column_name AS column_name,\n" +
                "    cc.comments AS comment\n" +
                "FROM all_col_comments cc\n" +
                "WHERE cc.owner = @schema\n" +
                "\n" +
                "ORDER BY object_type, object_name, column_name;\n";
    }
}

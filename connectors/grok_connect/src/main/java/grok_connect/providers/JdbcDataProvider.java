package grok_connect.providers;

import java.io.*;
import java.sql.*;
import java.util.*;
import java.math.*;
import java.text.*;
import java.util.regex.*;
import org.apache.commons.io.IOUtils;
import serialization.*;
import grok_connect.utils.*;
import grok_connect.table_query.*;
import grok_connect.connectors_info.*;
import serialization.Types;


public abstract class JdbcDataProvider extends DataProvider {
    public abstract Connection getConnection(DataConnection dataConnection)
            throws ClassNotFoundException, SQLException;

    public boolean isParametrized() {
        return true;
    }

    public String getConnectionString(DataConnection conn) {
        return conn.hasCustomConnectionString()
                ? (String)conn.parameters.get(DbCredentials.CONNECTION_STRING)
                : getConnectionStringImpl(conn);
    }

    public String getConnectionStringImpl(DataConnection conn) {
        return conn.connectionString;
    }

    // "CONN_AVAILABLE" or exception text
    public String testConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        try {
            Connection sqlConnection = getConnection(conn);
            if (sqlConnection.isClosed())
                return "Connection is not available";
            else
                return DataProvider.CONN_AVAILABLE;
        } catch (Throwable ex) {
            StringWriter errors = new StringWriter();
            errors.write("ERROR:\n" + ex.toString() + "\n\nSTACK TRACE:\n");
            ex.printStackTrace(new PrintWriter(errors));
            System.out.println(errors.toString());
            return errors.toString();
        }
    }

    public DataFrame getSchemas(DataConnection connection)
            throws ClassNotFoundException, SQLException, ParseException, IOException {
        FuncCall queryRun = new FuncCall();
        queryRun.func = new DataQuery();
        queryRun.func.query = getSchemasSql(connection.getDb());
        queryRun.func.connection = connection;

        return execute(queryRun);
    }

    public DataFrame getSchema(DataConnection connection, String schema, String table)
            throws ClassNotFoundException, SQLException, ParseException, IOException {
        FuncCall queryRun = new FuncCall();
        queryRun.func = new DataQuery();
        queryRun.func.query = getSchemaSql(connection.getDb(), schema, table);
        queryRun.func.connection = connection;

        return execute(queryRun);
    }

    public String getSchemasSql(String db) {
        throw new UnsupportedOperationException();
    }

    public String getSchemaSql(String db, String schema, String table) {
        throw new UnsupportedOperationException();
    }

    @SuppressWarnings("unchecked")
    public DataFrame execute(FuncCall queryRun)
            throws ClassNotFoundException, SQLException, ParseException, IOException {

        ResultSet resultSet;
        Pattern pattern = Pattern.compile("(?m)@(\\w+)");
        int count = (queryRun.options != null && queryRun.options.containsKey(DataProvider.QUERY_COUNT))
                ? ((Double)queryRun.options.get(DataProvider.QUERY_COUNT)).intValue() : 0;
        int memoryLimit = (queryRun.options != null && queryRun.options.containsKey(DataProvider.QUERY_MEMORY_LIMIT_MB))
                ? ((Double)queryRun.options.get(DataProvider.QUERY_MEMORY_LIMIT_MB)).intValue() : 0;
        int timeout = (queryRun.options != null && queryRun.options.containsKey(DataProvider.QUERY_TIMEOUT_SEC))
                ? ((Double)queryRun.options.get(DataProvider.QUERY_TIMEOUT_SEC)).intValue() : 300;

        // Remove header lines
        DataQuery dataQuery = queryRun.func;
        String query = dataQuery.query;
        String commentStart = DataProvider.getByName(dataQuery.connection.dataSource).descriptor.commentStart;
        query = query.replaceAll("(?m)^" + commentStart + ".*\\n", "");

        Connection connection = getConnection(dataQuery.connection);
        if (dataQuery.inputParamsCount() > 0) {
            query = convertPatternParamsToQueryParams(queryRun, query);

            if (isParametrized()) {
                // Parametrized func
                List<String> names = new ArrayList<>();
                Matcher matcher = pattern.matcher(query);
                StringBuilder queryBuffer = new StringBuilder();
                int idx = 0;
                while (matcher.find()) {
                    String name = matcher.group(1);
                    if (dataQuery.existsParam(name)) {
                        queryBuffer.append(query, idx, matcher.start());
                        queryBuffer.append("?");
                        idx = matcher.end();
                        names.add(name);
                    }
                }
                queryBuffer.append(query, idx, query.length());
                query = queryBuffer.toString();
                PreparedStatement statement = connection.prepareStatement(query);
                for (int n = 0; n < names.size(); n++) {
                    FuncParam param = dataQuery.getParam(names.get(n));
                    if (param.propertyType.equals(Types.DATE_TIME)) {
                        Calendar calendar = javax.xml.bind.DatatypeConverter.parseDateTime((String)param.value);
                        statement.setTimestamp(n + 1, new Timestamp(calendar.getTime().getTime()));
                    } else
                        statement.setObject(n + 1, param.value);
                }
                statement.setQueryTimeout(timeout);
                System.out.println(query);
                resultSet = statement.executeQuery();
            } else {
                // Put parameters into func
                Matcher matcher = pattern.matcher(query);
                StringBuilder queryBuffer = new StringBuilder();
                int idx = 0;
                while (matcher.find()) {
                    String name = matcher.group().substring(1);
                    queryBuffer.append(query.substring(idx, matcher.start()));
                    for (FuncParam param: dataQuery.getInputParams()) {
                        if (param.name.equals(name)) {
                            switch (param.propertyType) {
                                case Types.DATE_TIME:
                                    queryBuffer.append(castParamValueToSqlDateTime(param));
                                    break;
                                case Types.BOOL:
                                    queryBuffer.append(((boolean) param.value) ? "1=1" : "1=0");
                                    break;
                                case Types.STRING:
                                    queryBuffer.append("'");
                                    queryBuffer.append(param.value.toString());
                                    queryBuffer.append("'");
                                    break;
                                default:
                                    queryBuffer.append(param.value.toString());
                            }
                            break;
                        }
                    }
                    idx = matcher.end();
                }
                queryBuffer.append(query.substring(idx, query.length()));
                query = queryBuffer.toString();

                Statement statement = connection.createStatement();
                statement.setQueryTimeout(timeout);
                System.out.println(query);
                resultSet = statement.executeQuery(query);
            }
        } else {
            // Query without parameters
            Statement statement = connection.createStatement();
            statement.setQueryTimeout(timeout);
            System.out.println(query);
            resultSet = statement.executeQuery(query);
        }

        ResultSetMetaData resultSetMetaData = resultSet.getMetaData();

        int columnCount = resultSetMetaData.getColumnCount();
        List<Column> columns = new ArrayList<>(columnCount);
        List<Boolean> supportedType = new ArrayList<>(columnCount);
        List<Boolean> initColumn = new ArrayList<>(columnCount);
        for (int c = 1; c < columnCount + 1; c++) {
            Column column;

            int type = resultSetMetaData.getColumnType(c);
            String typeName = resultSetMetaData.getColumnTypeName(c);
            supportedType.add(c - 1, true);
            initColumn.add(c - 1, true);

            int precision = resultSetMetaData.getPrecision(c);
            int scale = resultSetMetaData.getScale(c);

            if (isInteger(type, typeName, precision, scale))
                column = new IntColumn();
            else if (isFloat(type, typeName) || isDecimal(type, typeName))
                column = new FloatColumn();
            else if (isBoolean(type, typeName))
                column = new BoolColumn();
            else if (isString(type, typeName) ||
                    typeName.equalsIgnoreCase("uuid") ||
                    typeName.equalsIgnoreCase("set"))
                column = new StringColumn();
            else if (isBigInt(type, typeName))
                column = new BigIntColumn();
            else if (isTime(type, typeName))
                column = new DateTimeColumn();
            else {
                column = new StringColumn();
                supportedType.set(c - 1, false);
                initColumn.set(c - 1, false);
            }

            column.name = resultSetMetaData.getColumnLabel(c);
            columns.add(c - 1, column);
        }

        int rowCount = 0;
        while (resultSet.next()) {
            rowCount++;

            for (int c = 1; c < columnCount + 1; c++) {
                Object value = resultSet.getObject(c);

                int type = resultSetMetaData.getColumnType(c);
                String typeName = resultSetMetaData.getColumnTypeName(c);
                int precision = resultSetMetaData.getPrecision(c);
                int scale = resultSetMetaData.getScale(c);

                if (supportedType.get(c - 1)) {
                    String colType = columns.get(c - 1).getType();
                    if (isInteger(type, typeName, precision, scale) || isBoolean(type, typeName) ||
                            colType.equals(Types.INT) || colType.equals(Types.BOOL))
                        if (value instanceof Short)
                            columns.get(c - 1).add(((Short)value).intValue());
                        else if (value instanceof Double)
                            columns.get(c - 1).add(((Double)value).intValue());
                        else if (value instanceof Float)
                            columns.get(c - 1).add(((Float)value).intValue());
                        else if (value instanceof BigDecimal)
                            columns.get(c - 1).add(((BigDecimal)value).intValue());
                        else
                            columns.get(c - 1).add(value);
                    else if (isString(type, typeName)) {
                        if (type == java.sql.Types.CLOB) {
                            Reader reader = ((Clob)value).getCharacterStream();
                            StringWriter writer = new StringWriter();
                            IOUtils.copy(reader, writer);
                            columns.get(c - 1).add(writer.toString());
                        } else
                            columns.get(c - 1).add(value);
                    } else if (isDecimal(type, typeName))
                        columns.get(c - 1).add((value == null) ? null : ((BigDecimal)value).floatValue());
                    else if (isFloat(type, typeName) || (colType.equals(Types.FLOAT)))
                        if (value instanceof Double)
                            columns.get(c - 1).add(new Float((Double)value));
                        else
                            columns.get(c - 1).add(value);
                    else if (isBigInt(type, typeName) ||
                            typeName.equalsIgnoreCase("uuid") ||
                            typeName.equalsIgnoreCase("set") ||
                            colType.equals(Types.STRING))
                        columns.get(c - 1).add((value != null) ? value.toString() : "");
                    else if (isTime(type, typeName)) {
                        java.util.Date time = (value instanceof java.sql.Timestamp)
                                ? java.util.Date.from(((java.sql.Timestamp)value).toInstant())
                                : ((java.util.Date) value);
                        columns.get(c - 1).add((time == null) ? null : time.getTime() * 1000.0);
                    }
                } else {
                    Column column = columns.get(c - 1);
                    if (!initColumn.get(c - 1) && value != null) {
                        if ((value instanceof Byte) || (value instanceof Short) || (value instanceof Integer)) {
                            column = new IntColumn();
                            column.addAll(new Integer[rowCount - 1]);
                        } else if ((value instanceof Float) || (value instanceof Double)) {
                            column = new FloatColumn();
                            column.addAll(new Float[rowCount - 1]);
                        } else if ((value instanceof Boolean)) {
                            column = new BoolColumn();
                            column.addAll(new Boolean[rowCount - 1]);
                        }
                        column.name = resultSetMetaData.getColumnLabel(c);
                        columns.set(c - 1, column);
                        initColumn.set(c - 1, true);

                        System.out.printf("Data type '%s' is not supported yet. Write as '%s'.\n",
                                resultSetMetaData.getColumnTypeName(c), column.getType());
                    }

                    if (value instanceof Double)
                        value = new Float((Double)value);
                    else if (!(value instanceof Byte) && !(value instanceof Short) &&
                            !(value instanceof Integer) && !(value instanceof Boolean))
                        value = (value != null) ? value.toString() : null;

                    column.add(value);
                }
            }

            if (rowCount % 1000 == 0) {
                int size = 0;
                for (Column column : columns)
                    size += column.memoryInBytes();
                size = ((count > 0) ? (int)((long)count * size / rowCount) : size) / 1000000;
                if (memoryLimit > 0 && size > memoryLimit)
                    throw new SQLException("Too large query result: " +
                            String.valueOf(size) + " > " + String.valueOf(memoryLimit) + " MB");
            }
        }

        connection.close();

        DataFrame dataFrame = new DataFrame();
        dataFrame.addColumns(columns);

        return dataFrame;
    }

    private static String paramToNamesString(FuncParam param, PatternMatcher matcher, String type,
                                             PatternMatcherResult result) {
        StringBuilder builder = new StringBuilder();
        for (int n = 0 ; n < matcher.values.size(); n++) {
            String name = param.name + String.valueOf(n);
            builder.append("@");
            builder.append(name);
            builder.append(",");
            result.params.add(new FuncParam(type, name, matcher.values.get(n)));
        }
        return builder.deleteCharAt(builder.length() - 1).toString();
    }

    public PatternMatcherResult numericPatternConverter(FuncParam param, PatternMatcher matcher) {
        PatternMatcherResult result = new PatternMatcherResult();
        String type = param.options.get("pattern");
        if (matcher.op.equals(PatternMatcher.NONE))
            result.query = "(1 = 1)";
        else if (matcher.op.equals(PatternMatcher.RANGE_NUM)) {
            String name0 = param.name + "R0";
            String name1 = param.name + "R1";
            result.query = "(" + matcher.colName + " >= @" + name0 + " AND " + matcher.colName + " <= @" + name1 + ")";
            result.params.add(new FuncParam(type, name0, matcher.values.get(0)));
            result.params.add(new FuncParam(type, name1, matcher.values.get(1)));
        } else if (matcher.op.equals(PatternMatcher.IN)) {
            String names = paramToNamesString(param, matcher, type, result);
            result.query = "(" + matcher.colName + " " + matcher.op + " (" + names + "))";
        } else {
            result.query = "(" + matcher.colName + " " + matcher.op + " @" + param.name + ")";
            result.params.add(new FuncParam(type, param.name, matcher.values.get(0)));
        }

        return result;
    }

    public PatternMatcherResult stringPatternConverter(FuncParam param, PatternMatcher matcher) {
        PatternMatcherResult result = new PatternMatcherResult();

        if (matcher.op.equals(PatternMatcher.NONE)) {
            result.query = "(1 = 1)";
            return result;
        }

        String type = "string";
        String _query = "(LOWER(" + matcher.colName + ") LIKE @" + param.name + ")";
        String value = ((String)matcher.values.get(0)).toLowerCase();

        if (matcher.op.equals(PatternMatcher.EQUALS)) {
            result.query = _query;
            result.params.add(new FuncParam(type, param.name, value));
        } else if (matcher.op.equals(PatternMatcher.CONTAINS)) {
            result.query = _query;
            result.params.add(new FuncParam(type, param.name, "%" + value + "%"));
        } else if (matcher.op.equals(PatternMatcher.STARTS_WITH)) {
            result.query = _query;
            result.params.add(new FuncParam(type, param.name, value + "%"));
        } else if (matcher.op.equals(PatternMatcher.ENDS_WITH)) {
            result.query = _query;
            result.params.add(new FuncParam(type, param.name, "%" + value));
        } else if (matcher.op.equals(PatternMatcher.REGEXP)) {
            result.query = "(1 = 1)"; // TODO Implement regexp
        } else if (matcher.op.equals(PatternMatcher.IN)) {
            String names = paramToNamesString(param, matcher, type, result);
            result.query = "(" + matcher.colName + " " + matcher.op + " (" + names + "))";
        } else {
            result.query = "(1 = 1)";
        }

        return result;
    }

    public PatternMatcherResult dateTimePatternConverter(FuncParam param, PatternMatcher matcher) {
        PatternMatcherResult result = new PatternMatcherResult();

        if (matcher.op.equals(PatternMatcher.NONE)) {
            result.query = "(1 = 1)";
        } else if (matcher.op.equals(PatternMatcher.EQUALS)) {
            result.query = "(" + matcher.colName + " = @" + param.name + ")";
            result.params.add(new FuncParam("datetime", param.name, matcher.values.get(0)));
        } else if (matcher.op.equals(PatternMatcher.BEFORE) || matcher.op.equals(PatternMatcher.AFTER)) {
            result.query = "(" + matcher.colName + PatternMatcher.cmp(matcher.op, matcher.include1) + "@" + param.name + ")";
            result.params.add(new FuncParam("datetime", param.name, matcher.values.get(0)));
        } else if (matcher.op.equals(PatternMatcher.RANGE_DATE_TIME)) {
            String name0 = param.name + "R0";
            String name1 = param.name + "R1";
            result.query = "(" + matcher.colName + PatternMatcher.cmp(PatternMatcher.AFTER, matcher.include1) + "@" + name0 + " AND " +
                    matcher.colName + PatternMatcher.cmp(PatternMatcher.BEFORE, matcher.include2) + "@" + name1 + ")";
            result.params.add(new FuncParam("datetime", name0, matcher.values.get(0)));
            result.params.add(new FuncParam("datetime", name1, matcher.values.get(1)));
        } else {
            result.query = "(1 = 1)";
        }

        return result;
    }

    private String aggrToSql(GroupAggregation aggr) {
        AggrFunctionInfo funcInfo = null;
        for (AggrFunctionInfo info: descriptor.aggregations) {
            if (info.functionName.equals(aggr.aggType)) {
                funcInfo = info;
                break;
            }
        }
        if (funcInfo != null) {
            String sql = funcInfo.dbFunctionName.replaceAll("#", aggr.colName);
            return sql + " as \"" + sql + "\"";
        } else
            return null;
    }

    public String addBrackets(String name) {
        String brackets = descriptor.nameBrackets;
        return (name.contains(" ") && !name.startsWith(brackets.substring(0, 1)))
                ? brackets.substring(0, 1) + name + brackets.substring(brackets.length() - 1, brackets.length())
                : name;
    }

    public String limitToSql(String query, Integer limit) {
        return query + "limit " + limit.toString();
    }

    private String patternToSql(FieldPredicate condition) {
        return condition.matcher.toSql(condition.dataType, condition.field);
    }

    public String queryTableSql(DataConnection conn, TableQuery query) {
        return query.toSql(this::aggrToSql, this::patternToSql, this::limitToSql, this::addBrackets,
                descriptor.limitAtEnd);
    }

    public DataFrame queryTable(DataConnection conn, TableQuery query)
            throws ClassNotFoundException, SQLException, ParseException, IOException {
        FuncCall queryRun = new FuncCall();
        queryRun.func = new DataQuery();
        String sql = queryTableSql(conn, query);
        System.out.println(sql);
        if (sql == null)
            return new DataFrame();
        queryRun.func.query = sql;
        queryRun.func.connection = conn;
        return execute(queryRun);
    }

    public String castParamValueToSqlDateTime(FuncParam param) {
        return "datetime('" + param.value.toString() + "')";
    }

    // TODO Convert following code into "List.contains() style
    private static boolean isInteger(int type, String typeName, int precision, int scale) {
        return (type == java.sql.Types.INTEGER) || (type == java.sql.Types.TINYINT) || (type == java.sql.Types.SMALLINT) ||
                typeName.equalsIgnoreCase("int4") ||
                typeName.equalsIgnoreCase("int2") ||
                typeName.equalsIgnoreCase("int") ||
                typeName.equalsIgnoreCase("serial2") ||
                typeName.equalsIgnoreCase("serial4") ||
                ((precision < 33) && (scale == 0) && (isFloat(type, typeName) || isDecimal(type, typeName)));
        // TODO Investigate precision value for current case
    }

    private static boolean isTime(int type, String typeName) {
        return (type == java.sql.Types.DATE) || (type == java.sql.Types.TIME) || (type == java.sql.Types.TIMESTAMP) ||
                typeName.equalsIgnoreCase("timetz") ||
                typeName.equalsIgnoreCase("timestamptz");
    }

    private static boolean isBigInt(int type, String typeName) {
        return (type == java.sql.Types.BIGINT) || typeName.equalsIgnoreCase("int8") ||
                typeName.equalsIgnoreCase("serial8");
    }

    private static boolean isFloat(int type, String typeName) {
        return (type == java.sql.Types.FLOAT) || (type == java.sql.Types.DOUBLE) || (type == java.sql.Types.REAL) ||
                typeName.equalsIgnoreCase("float8") ||
                typeName.equalsIgnoreCase("float4") ||
                typeName.equalsIgnoreCase("money");
    }

    private static boolean isDecimal(int type, String typeName) {
        return (type == java.sql.Types.DECIMAL) || (type == java.sql.Types.NUMERIC) ||
                typeName.equalsIgnoreCase("decimal");
    }

    private static boolean isBoolean(int type, String typeName) {
        return (type == java.sql.Types.BOOLEAN) || (type == java.sql.Types.BIT) ||
                typeName.equalsIgnoreCase("bool");
    }

    private static boolean isString(int type, String typeName) {
        return ((type == java.sql.Types.VARCHAR)|| (type == java.sql.Types.CHAR) ||
                (type == java.sql.Types.LONGVARCHAR) || (type == java.sql.Types.CLOB) ||
                typeName.equalsIgnoreCase("varchar") ||
                typeName.equalsIgnoreCase("nvarchar") ||
                typeName.equalsIgnoreCase("nchar") ||
                typeName.equalsIgnoreCase("ntext")) &&
                !typeName.equalsIgnoreCase("uuid") &&
                !typeName.equalsIgnoreCase("set");
    }

    public static java.util.Properties defaultConnectionProperties(DataConnection conn) {
        java.util.Properties properties = new java.util.Properties();
        properties.setProperty("user", conn.credentials.getLogin());
        properties.setProperty("password", conn.credentials.getPassword());
        return properties;
    }
}

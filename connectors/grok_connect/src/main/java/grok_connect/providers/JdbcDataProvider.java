package grok_connect.providers;

import java.sql.Array;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;
import java.sql.Statement;
import java.sql.Timestamp;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataProvider;
import grok_connect.connectors_info.DataQuery;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncCall;
import grok_connect.connectors_info.FuncParam;
import grok_connect.log.EventType;
import grok_connect.resultset.DefaultResultSetManager;
import grok_connect.resultset.ResultSetManager;
import grok_connect.table_mutation.AddKey;
import grok_connect.table_mutation.AlterTable;
import grok_connect.table_mutation.BatchInsertBulkLoader;
import grok_connect.table_mutation.BulkLoader;
import grok_connect.table_mutation.ColumnSpec;
import grok_connect.table_mutation.CreateIndex;
import grok_connect.table_mutation.CreateTable;
import grok_connect.table_mutation.DeleteRows;
import grok_connect.table_mutation.DropIndex;
import grok_connect.table_mutation.DropTable;
import grok_connect.table_mutation.IndexSpec;
import grok_connect.table_mutation.InsertRows;
import grok_connect.table_mutation.MutationValidationException;
import grok_connect.table_mutation.TableMutation;
import grok_connect.table_mutation.TruncateTable;
import grok_connect.table_mutation.UpdateRows;
import grok_connect.table_mutation.UpsertRows;
import grok_connect.table_query.AggrFunctionInfo;
import grok_connect.table_query.FieldPredicate;
import grok_connect.table_query.GroupAggregation;
import grok_connect.table_query.PredicateCompiler;
import grok_connect.table_query.SqlNames;
import grok_connect.table_query.TableQuery;
import grok_connect.utils.*;
import org.apache.commons.lang.NotImplementedException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.DataFrame;
import serialization.StringColumn;
import serialization.Types;

public abstract class JdbcDataProvider extends DataProvider {
    public static Pattern UUID_REGEX = Pattern.compile("^[\\da-fA-F]{8}-[\\da-fA-F]{4}-[\\da-fA-F]{4}-[\\da-fA-F]{4}-[\\da-fA-F]{12}$");
    protected Logger logger = LoggerFactory.getLogger(this.getClass().getName());
    protected QueryMonitor queryMonitor = QueryMonitor.getInstance();
    protected String driverClassName;

    public void prepareProvider() throws GrokConnectException {
    }

    public Connection getConnection(DataConnection conn) throws SQLException, GrokConnectException {
        prepareProvider();
        return ConnectionPool.getConnection(getConnectionString(conn), getProperties(conn), driverClassName);
    }

    public Properties getProperties(DataConnection conn) {
        java.util.Properties props = defaultConnectionProperties(conn);
        java.util.Properties jdbcProps = getJdbcProperties(conn);
        props.putAll(jdbcProps);
        return props;
    }

    @SuppressWarnings("unchecked")
    public Properties getJdbcProperties(DataConnection conn) {
        java.util.Properties props = new Properties();
        if (conn.parameters.containsKey("jdbcProperties") && this.descriptor.jdbcPropertiesTemplate != null) {
            Map<String, Object> propMap = (Map<String, Object>) conn.parameters.get("jdbcProperties");
            if (!propMap.isEmpty()) {
                for (Property p : this.descriptor.jdbcPropertiesTemplate) {
                    Object value = propMap.get(p.name);
                    if (value != null) {
                        String strValue = value.toString();
                        if (GrokConnectUtil.isEmpty(strValue))
                            continue;
                        if (p.propertyType.equals(Property.INT_TYPE))
                            strValue = String.valueOf(new Double(strValue).intValue());
                        props.setProperty(p.name, strValue);
                    }
                }
            }
        }
        return props;
    }

    public boolean autoInterpolation() {
        return true;
    }

    protected Integer getTimeout() {
        return null;
    }

    public String getConnectionString(DataConnection conn) {
        return conn.hasCustomConnectionString()
                ? (String)conn.parameters.get(DbCredentials.CONNECTION_STRING)
                : getConnectionStringImpl(conn);
    }

    public String getConnectionStringImpl(DataConnection conn) {
        return conn.connectionString;
    }

    public void testConnection(DataConnection conn) throws GrokConnectException {
        try (Connection ignored = getConnection(conn)) {
            // just open and close the connection
        } catch (SQLException e) {
            throw new GrokConnectException(e);
        }
    }

    public DataFrame getSchemas(DataConnection connection) throws QueryCancelledByUser, GrokConnectException {
        FuncCall queryRun = new FuncCall();
        queryRun.func = new DataQuery();
        queryRun.func.query = getSchemasSql(connection.getDb());
        queryRun.func.connection = connection;

        return execute(queryRun);
    }

    public DataFrame getSchema(DataConnection connection, String schema, String table, boolean includeKeyInfo) throws QueryCancelledByUser,
            GrokConnectException {
        FuncCall queryRun = new FuncCall();
        queryRun.func = new DataQuery();
        queryRun.func.query = getSchemaSql(connection.getDb(), schema, table, includeKeyInfo);
        queryRun.func.connection = connection;

        return execute(queryRun);
    }

    public DataFrame getCatalogs(DataConnection connection) throws GrokConnectException {
        try (Connection db = getConnection(connection)) {
            return readCatalogsFromMetadata(db);
        } catch (SQLException e) {
            throw new GrokConnectException(e);
        }
    }

    protected DataFrame readCatalogsFromMetadata(Connection connection) throws SQLException {
        DataFrame result = DataFrame.fromColumns(new StringColumn("catalog_name"));
        try (ResultSet rs = connection.getMetaData().getCatalogs()) {
            List<String> catalogs = new ArrayList<>();
            while (rs.next())
                catalogs.add(rs.getString("TABLE_CAT"));
            Collections.sort(catalogs);
            for (String catalog : catalogs)
                result.addRow(catalog);
        }
        return result;
    }

    public String getCommentsQuery(DataConnection connection) throws GrokConnectException {
        throw new GrokConnectException("Operation is not supported");
    }

    public DataFrame getComments(DataConnection connection, String schema) throws GrokConnectException, QueryCancelledByUser {
        FuncCall queryRun = new FuncCall();
        queryRun.func = new DataQuery();
        queryRun.func.query = getCommentsQuery(connection);
        queryRun.func.parameters.put("schema", schema);
        queryRun.func.params = new ArrayList<>();
        queryRun.func.params.add(new FuncParam("string", "schema", schema));
        queryRun.func.connection = connection;

        return execute(queryRun);
    }

    public List<String> getUniqueColumns(DataConnection conn, String schema, String table) throws GrokConnectException {
        try (Connection connection = getConnection(conn);
             ResultSet rs = connection.getMetaData().getIndexInfo(conn.getDb(), schema, table, true, false)) {
            List<String> uniqueColumns = new ArrayList<>();
            while(rs.next())
                uniqueColumns.add(rs.getString("COLUMN_NAME"));
            return uniqueColumns;

        } catch (SQLException e) {
            throw new GrokConnectException(e);
        }
    }

    public DataFrame getForeignKeys(DataConnection conn, String schema) throws GrokConnectException, QueryCancelledByUser {
        try (Connection connection = getConnection(conn);
             ResultSet rs = connection.getMetaData().getExportedKeys(null, schema, null)) {
            DataFrame result = DataFrame.fromColumns(new StringColumn("table_schema"),
                    new StringColumn("constraint_name"), new StringColumn("table_name"),
                    new StringColumn("column_name"), new StringColumn("foreign_table_name"), new StringColumn("foreign_column_name"));
            while(rs.next())
                result.addRow(rs.getString("FKTABLE_SCHEM"), rs.getString("FK_NAME"), rs.getString("FKTABLE_NAME"),
                        rs.getString("FKCOLUMN_NAME"), rs.getString("PKTABLE_NAME"), rs.getString("PKCOLUMN_NAME"));
            return result;

        } catch (SQLException e) {
            throw new GrokConnectException(e);
        }
    }

    public String getSchemasSql(String db) {
        throw new UnsupportedOperationException();
    }

    public String getSchemaSql(String db, String schema, String table, boolean includeKeyInfo) {
        throw new UnsupportedOperationException();
    }

    public ResultSet executeQuery(String query, FuncCall queryRun,
                                  Connection connection, int timeout, int fetchSize) throws SQLException {
        DataQuery dataQuery = queryRun.func;
        String mainCallId = (String) queryRun.aux.get("mainCallId");

        ResultSet resultSet;
        if (dataQuery.inputParamsCount() > 0) {
            logger.debug(EventType.QUERY_PARSE.getMarker(EventType.Stage.START), "Converting query pattern parameters...");
            query = convertPatternParamsToQueryParams(queryRun, query);
            logger.debug(EventType.QUERY_PARSE.getMarker(EventType.Stage.END), "Converted query pattern parameters");
            if (autoInterpolation()) {
                StringBuilder queryBuffer = new StringBuilder();
                logger.debug(EventType.QUERY_INTERPOLATION.getMarker(EventType.Stage.START), "Interpolating automatically SQL query parameters...");
                List<String> names = getParameterNames(query, dataQuery, queryBuffer);
                query = queryBuffer.toString();
                logger.debug(EventType.QUERY_INTERPOLATION.getMarker(EventType.Stage.END), "Interpolated SQL query parameters. Detected {} parameters", names.size());
                logger.info("Query before execution: {}", query);

                logger.debug("Creating PreparedStatement...");
                PreparedStatement statement = connection.prepareStatement(query);
                logger.debug("Created PreparedStatement");
                logger.debug(EventType.STATEMENT_PARAMETERS_REPLACEMENT.getMarker(EventType.Stage.START), "Replacing designated query parameters ? with actual values...");
                int i = 0;
                for (int n = 0; n < names.size(); n++) {
                    FuncParam param = dataQuery.getParam(names.get(n));
                    if (param.propertyType.equals(Types.DATE_TIME)) {
                        setDateTimeValue(param, statement, n + i + 1);
                    }
                    else if (param.propertyType.equals(Types.LIST) && param.propertySubType.equals(Types.STRING)) {
                        i = i + setArrayParamValue(statement, n + i + 1, param);
                    }
                    else if (param.propertyType.equals(Types.BIG_INT) && param.value != null)
                        statement.setLong(n + i + 1, Long.parseLong(param.value.toString()));
                    else if (param.propertyType.equals(Types.STRING) && param.value != null) {
                        String s = param.value.toString();
                        if (UUID_REGEX.matcher(s).matches())
                            setUuid(statement, n + i + 1, s);
                        else
                            statement.setString(n + i + 1, s);
                    }
                    else {
                        if (param.value == null) {
                            switch (param.propertyType) {
                                case Types.INT:
                                case Types.FLOAT:
                                    statement.setNull(n + i + 1, java.sql.Types.NUMERIC);
                                    break;
                                case Types.BIG_INT:
                                    statement.setNull(n + i + 1, java.sql.Types.BIGINT);
                                case Types.BOOL:
                                    statement.setNull(n + i + 1, java.sql.Types.BOOLEAN);
                                default:
                                    statement.setNull(n + i + 1, java.sql.Types.VARCHAR);
                                    break;
                            }
                        }
                        else
                            statement.setObject(n + i + 1, param.value);
                    }
                }
                logger.debug(EventType.STATEMENT_PARAMETERS_REPLACEMENT.getMarker(EventType.Stage.END), "Replaced designated query parameters");
                resultSet = executeStatement(statement, timeout, mainCallId, fetchSize);
            } else {
                logger.debug(EventType.QUERY_INTERPOLATION.getMarker(EventType.Stage.START), "Interpolating manually SQL query parameters...");
                query = manualQueryInterpolation(query, dataQuery);
                logger.debug(EventType.QUERY_INTERPOLATION.getMarker(EventType.Stage.END), "Interpolated SQL query parameters");
                logger.info("Query before execution: {}", query);
                resultSet = executeStatement(connection.prepareStatement(query), timeout, mainCallId, fetchSize);
            }
        } else {
            logger.info("Query before execution: {}", query);
            resultSet = executeStatement(connection.prepareStatement(query), timeout, mainCallId, fetchSize);
        }

        return resultSet;
    }

    private ResultSet executeStatement(PreparedStatement statement, int timeout, String mainCallId,
                                       int fetchSize) throws SQLException {
        queryMonitor.addNewStatement(mainCallId, statement);
        setQueryTimeOut(statement, timeout);
        logger.debug(EventType.STATEMENT_EXECUTION.getMarker(EventType.Stage.START), "Executing Statement...");
        statement.setFetchSize(fetchSize);
        ResultSet resultSet = executeStatement(statement);
        logger.info(EventType.STATEMENT_EXECUTION.getMarker(EventType.Stage.END), "Executed Statement");
        queryMonitor.removeStatement(mainCallId);
        return resultSet;
    }

    protected ResultSet executeStatement(PreparedStatement statement) throws SQLException {
        return statement.execute() ? statement.getResultSet() : null;
    }

    private void setQueryTimeOut(Statement statement, int timeout) {
        try {
            statement.setQueryTimeout(timeout);
        } catch (SQLException exception) {
            logger.debug("setQueryTimeout is not supported for {}", descriptor.type);
        }
    }

    public void setDateTimeValue(FuncParam funcParam, PreparedStatement statement, int parameterIndex) throws SQLException {
        if (funcParam.value == null) {
            statement.setNull(parameterIndex, java.sql.Types.TIMESTAMP);
            return;
        }
        Calendar calendar = javax.xml.bind.DatatypeConverter.parseDateTime((String)funcParam.value);
        calendar.setTimeZone(TimeZone.getTimeZone("UTC"));
        Timestamp ts = new Timestamp(calendar.getTime().getTime());
        statement.setTimestamp(parameterIndex, ts, Calendar.getInstance(TimeZone.getTimeZone("UTC")));
    }

    protected String manualQueryInterpolation(String query, DataQuery dataQuery) {
        Pattern pattern = Pattern.compile("(?m)@(\\w+)");
        // Put parameters into func
        Matcher matcher = pattern.matcher(query);
        StringBuilder queryBuffer = new StringBuilder();
        int idx = 0;
        while (matcher.find()) {
            String name = matcher.group().substring(1);
            queryBuffer.append(query, idx, matcher.start());
            interpolateParameters(queryBuffer, dataQuery, name);
            idx = matcher.end();
        }
        queryBuffer.append(query.substring(idx));
        query = queryBuffer.toString();
        return query;
    }

    protected void interpolateParameters(StringBuilder queryBuffer, DataQuery dataQuery, String paramName) {
        for (FuncParam param: dataQuery.getInputParams()) {
            if (param.name.equals(paramName)) {
                switch (param.propertyType) {
                    case Types.DATE_TIME:
                        queryBuffer.append(castParamValueToSqlDateTime(param));
                        return;
                    case Types.BOOL:
                        queryBuffer.append(interpolateBool(param));
                        return;
                    case Types.STRING: //todo: support escaping
                        queryBuffer.append(interpolateString(param));
                        return;
                    case Types.LIST: //todo: extract submethod
                        if (param.propertySubType.equals(Types.STRING)) {
                            @SuppressWarnings(value = "unchecked")
                            ArrayList<String> value = ((ArrayList<String>) param.value);
                            for (int i = 0; i < value.size(); i++) {
                                queryBuffer.append(String.format("'%s'", value.get(i)));
                                if (i < value.size() - 1)
                                    queryBuffer.append(",");
                            }
                            return;
                            //todo: implement other types
                        } else {
                            throw new NotImplementedException("Non-string lists are not implemented for manual param interpolation providers");
                        }
                    default:
                        queryBuffer.append(param.value.toString());
                }
                return;
            }
        }
        queryBuffer
                .append("@")
                .append(paramName); // there are no such FuncParam, so it means that it is not a param
    }

    protected String interpolateString(FuncParam param) {
        return String.format("'%s'", param.value.toString());
    }

    protected String interpolateBool(FuncParam param) {
        return ((boolean) param.value) ? "1=1" : "1=0";
    }

    protected int setArrayParamValue(PreparedStatement statement, int n, FuncParam param) throws SQLException {
        if (param.value == null)
            statement.setNull(n, java.sql.Types.ARRAY);
        else {
            @SuppressWarnings (value="unchecked")
            ArrayList<String> values = (ArrayList<String>) param.value;
            Array array = statement.getConnection().createArrayOf("VARCHAR", values.toArray());
            statement.setArray(n, array);
        }
        return 0;
    }

    protected void setUuid(PreparedStatement statement, int n, String value) throws SQLException {
        statement.setString(n, value);
    }

    public List<String> getParameterNames(String query, DataQuery dataQuery, StringBuilder queryBuffer) {
        List<String> names = new ArrayList<>();
        String regexComment = String.format("(?m)^(?<!['\\\"])%s.*(?!['\\\"])$", descriptor.commentStart);
        query = query
                .replaceAll(regexComment, "")
                .trim();
        Pattern pattern = Pattern.compile("(?m)@(\\w+)");
        Matcher matcher = pattern.matcher(query);
        int idx = 0;
        while (matcher.find()) {
            String name = matcher.group(1);
            if (dataQuery.existsParam(name)) {
                queryBuffer.append(query, idx, matcher.start());
                appendQueryParam(dataQuery, name, queryBuffer);
                idx = matcher.end();
                names.add(name);
            }
        }
        queryBuffer.append(query, idx, query.length());
        return names;
    }

    protected void appendQueryParam(DataQuery dataQuery, String paramName, StringBuilder queryBuffer) {
        queryBuffer.append("?");
    }

    public void configureAutoCommit(Connection connection) throws SQLException {
        boolean supportsTransactions = connection.getMetaData().supportsTransactions();
        logger.trace("Provider {} transactions", supportsTransactions ? "supports" : "doesn't support");
        if (supportsTransactions)
            connection.setAutoCommit(false);
    }

    public ResultSet getResultSet(FuncCall queryRun, Connection connection, int fetchSize) throws QueryCancelledByUser, SQLException {
        Integer providerTimeout = getTimeout();
        int timeout = providerTimeout != null ? providerTimeout : (queryRun.options != null && queryRun.options.containsKey(DataProvider.QUERY_TIMEOUT_SEC))
                ? ((Double)queryRun.options.get(DataProvider.QUERY_TIMEOUT_SEC)).intValue() : 300;
        configureAutoCommit(connection);
        try {
            // Remove header lines
            DataQuery dataQuery = queryRun.func;
            String query = dataQuery.query;
            String commentStart = descriptor.commentStart;

            ResultSet resultSet = null;

            if (!(queryRun.func.options != null
                    && queryRun.func.options.containsKey("batchMode")
                    && queryRun.func.options.get("batchMode").equals("true"))) {
                query = query.replaceAll("(?m)^\\s*" + commentStart + ".*\\n", "");
                resultSet = executeQuery(query, queryRun, connection, timeout, fetchSize);
            }
            else {
                logger.debug("Executing batch mode...");
                String[] queries = query.replaceAll("\r\n", "\n").split(String.format("\n%sbatch\n|\n--batch\n", commentStart));
                for (String currentQuery : queries)
                    resultSet = executeQuery(currentQuery, queryRun, connection, timeout, fetchSize);
                logger.debug("Executed batch mode");
            }

            return resultSet;
        } catch (SQLException e) {
            if (queryMonitor.checkCancelledId((String) queryRun.aux.get("mainCallId")))
                throw new QueryCancelledByUser();
            else throw e;
        }
    }
    public DataFrame getResultSetSubDf(FuncCall queryRun, ResultSet resultSet, ResultSetManager resultSetManager,
                                       int maxIterations, int columnCount, int operationNumber,
                                       boolean dryRun) throws SQLException, QueryCancelledByUser {
        DataFrame dataFrame = new DataFrame();
        EventType resultSetProcessingEventType = dryRun ? EventType.RESULT_SET_PROCESSING_WITHOUT_DATAFRAME_FILL
                : EventType.RESULT_SET_PROCESSING_WITH_DATAFRAME_FILL;
        logger.debug(resultSetProcessingEventType.getMarker(operationNumber, EventType.Stage.START),
                "Filling columns of DataFrame with id {}...", operationNumber);
        if (resultSet.next()) {
            int rowCount = 0;
            do {
                rowCount++;
                for (int c = 1; c < columnCount + 1; c++) {
                    Object value = getObjectFromResultSet(resultSet, c);

                    if (dryRun) continue;
                    resultSetManager.processValue(value, c);

                    if (queryMonitor.checkCancelledIdResultSet(queryRun.id)) {
                        logger.info("Query was canceled");
                        queryMonitor.removeResultSet(queryRun.id);
                        throw new QueryCancelledByUser();
                    }
                }
            } while ((maxIterations < 0 || rowCount < maxIterations) && resultSet.next());

            logger.debug(resultSetProcessingEventType.getMarker(operationNumber, EventType.Stage.END),
                    "Filled columns with {} rows of DataFrame with id {}", rowCount, operationNumber);
        }
        else
            logger.debug(resultSetProcessingEventType.getMarker(operationNumber, EventType.Stage.END),
                    "Result set is empty");
        dataFrame.addColumns(resultSetManager.getProcessedColumns());
        return dataFrame;
    }

    public DataFrame execute(FuncCall queryRun) throws QueryCancelledByUser, GrokConnectException {
        if (queryRun.func instanceof TableMutation)
            throw new GrokConnectException("Table mutations must be executed via POST /mutate, not the query path");
        Connection connection = null;
        ResultSet resultSet = null;
        try {
            connection = getConnection(queryRun.func.connection);
            resultSet = getResultSet(queryRun, connection, 100);
            DataFrame dataFrame = new DataFrame();
            if (resultSet != null) {
                ResultSetManager resultSetManager = getResultSetManager();
                ResultSetMetaData metaData = resultSet.getMetaData();
                resultSetManager.init(metaData, 100);
                dataFrame = getResultSetSubDf(queryRun, resultSet, resultSetManager, -1, metaData.getColumnCount(), 1, false);
            }
            if (!connection.getAutoCommit())
                connection.commit();
            return dataFrame;
        } catch (SQLException e) {
            rollbackQuietly(connection);
            if (queryMonitor.checkCancelledId((String) queryRun.aux.get("mainCallId")))
                throw new QueryCancelledByUser();
            else throw new GrokConnectException(e);
        } catch (QueryCancelledByUser | RuntimeException e) {
            rollbackQuietly(connection);
            throw e;
        } finally {
            if (resultSet != null)
                try {
                    resultSet.close();
                } catch (SQLException e) {
                    logger.warn("Failed to close ResultSet", e);
                }
            if (connection != null)
                try {
                    connection.close();
                } catch (SQLException e) {
                    logger.warn("Failed to close connection", e);
                }
        }
    }

    public void rollbackQuietly(Connection connection) {
        if (connection == null)
            return;
        try {
            if (!connection.getAutoCommit())
                connection.rollback();
        } catch (SQLException e) {
            logger.warn("Failed to rollback transaction", e);
        }
    }

    protected Object getObjectFromResultSet(ResultSet resultSet, int c) {
        try {
            return resultSet.getObject(c);
        }catch (SQLException e) {
            throw new RuntimeException("Something went wrong when getting object from result set", e);
        }
    }

    protected static String paramToNamesString(String paramName, PatternMatcher matcher, String type,
                                               PatternMatcherResult result) {
        StringBuilder builder = new StringBuilder();
        for (int n = 0 ; n < matcher.values.size(); n++) {
            String name = paramName + n;
            builder.append("@");
            builder.append(name);
            builder.append(",");
            result.params.add(new FuncParam(type, name, matcher.values.get(n)));
        }
        return builder.deleteCharAt(builder.length() - 1).toString();
    }

    @Override
    public PatternMatcherResult numericPatternConverter(String paramName, String typeName, PatternMatcher matcher) {
        PatternMatcherResult result = new PatternMatcherResult();
        switch (matcher.op) {
            case PatternMatcher.NONE:
                result.query = "(1 = 1)";
                break;
            case PatternMatcher.RANGE_NUM:
                String name0 = paramName + "R0";
                String name1 = paramName + "R1";
                result.query = "(" + matcher.colName + " >= @" + name0 + " AND " + matcher.colName + " <= @" + name1 + ")";
                result.params.add(new FuncParam(typeName, name0, matcher.values.get(0)));
                result.params.add(new FuncParam(typeName, name1, matcher.values.get(1)));
                break;
            case PatternMatcher.IN:
            case PatternMatcher.NOT_IN:
                String names = paramToNamesString(paramName, matcher, typeName, result);
                result.query = getInQuery(matcher, names);
                break;
            case PatternMatcher.IS_NULL:
            case PatternMatcher.IS_NOT_NULL:
                result.query = String.format("(%s %s)", matcher.colName, matcher.op);
                break;
            default:
                result.query = "(" + matcher.colName + " " + matcher.op + " @" + paramName + ")";
                result.params.add(new FuncParam(typeName, paramName, matcher.values.get(0)));
                break;
        }
        return result;
    }

    protected String getInQuery(PatternMatcher matcher, String names) {
        return String.format("(%s %s (%s))", matcher.colName, matcher.op, names);
    }

    @Override
    public PatternMatcherResult boolPatternConverter(String paramName, PatternMatcher matcher) {
        PatternMatcherResult result = new PatternMatcherResult();
        if (matcher.op.equals(PatternMatcher.EQUALS)) {
            result.query = "(" + matcher.colName + " = @" + paramName + ")";
            result.params.add(new FuncParam(Types.BOOL, paramName, matcher.values.stream().findFirst().orElse("true")));
        }
        else
            result.query = "(1 = 1)";
        return result;
    }

    @Override
    public PatternMatcherResult bigIntPatternConverter(String paramName, PatternMatcher matcher) {
        PatternMatcherResult result = new PatternMatcherResult();
        if (matcher.op.equals(PatternMatcher.EQUALS)) {
            result.query = "(" + matcher.colName + " = @" + paramName + ")";
            result.params.add(new FuncParam(Types.BIG_INT, paramName, matcher.values.stream().findFirst().orElse(null)));
        }
        else if (matcher.op.equals(PatternMatcher.IN) || matcher.op.equals(PatternMatcher.NOT_IN)) {
            String names = paramToNamesString(paramName, matcher, "bigint", result);
            result.query = getInQuery(matcher, names);
        }
        else
            result.query = "(1 = 1)";
        return result;
    }

    @Override
    public PatternMatcherResult stringPatternConverter(String paramName, PatternMatcher matcher) {
        PatternMatcherResult result = new PatternMatcherResult();

        if (matcher.op.equals(PatternMatcher.NONE)) {
            result.query = "(1 = 1)";
            return result;
        }

        String type = "string";
        String _query = "(LOWER(" + matcher.colName + ") LIKE @" + paramName + ")";
        List<Object> values = matcher.values;
        String value = null;
        if (values.size() > 0)
            value = ((String) values.get(0)).toLowerCase();

        switch (matcher.op) {
            case PatternMatcher.EQUALS:
                result.query = _query;
                result.params.add(new FuncParam(type, paramName, value));
                break;
            case PatternMatcher.CONTAINS:
                result.query = _query;
                result.params.add(new FuncParam(type, paramName, "%" + value + "%"));
                break;
            case PatternMatcher.STARTS_WITH:
                result.query = _query;
                result.params.add(new FuncParam(type, paramName, value + "%"));
                break;
            case PatternMatcher.ENDS_WITH:
                result.query = _query;
                result.params.add(new FuncParam(type, paramName, "%" + value));
                break;
            case PatternMatcher.REGEXP:
                result.query = getRegexQuery(matcher.colName, value);
                result.params.add(new FuncParam(type, paramName, value));
                break;
            case PatternMatcher.IN:
            case PatternMatcher.NOT_IN:
                String names = paramToNamesString(paramName, matcher, type, result);
                result.query = getInQuery(matcher, names);
                break;
            case PatternMatcher.IS_NULL:
            case PatternMatcher.IS_NOT_NULL:
                result.query = String.format("(%s %s)", matcher.colName, matcher.op);
                break;
            default:
                result.query = "(1 = 1)";
                break;
        }

        return result;
    }

    protected String getRegexQuery(String columnName, String regexExpression) {
        throw new UnsupportedOperationException("REGEXP is not supported for this provider");
    }

    @Override
    public PatternMatcherResult dateTimePatternConverter(String paramName, PatternMatcher matcher) {
        PatternMatcherResult result = new PatternMatcherResult();

        switch (matcher.op) {
            case PatternMatcher.EQUALS:
                result.query = "(" + matcher.colName + " = @" + paramName + ")";
                result.params.add(new FuncParam("datetime", paramName, matcher.values.get(0)));
                break;
            case PatternMatcher.BEFORE:
            case PatternMatcher.AFTER:
                result.query = "(" + matcher.colName + PatternMatcher.cmp(matcher.op, matcher.include1) + "@" + paramName + ")";
                result.params.add(new FuncParam("datetime", paramName, matcher.values.get(0)));
                break;
            case PatternMatcher.RANGE_DATE_TIME:
                String name0 = paramName + "R0";
                String name1 = paramName + "R1";
                result.query = "(" + matcher.colName + PatternMatcher.cmp(PatternMatcher.AFTER, matcher.include1) + "@" + name0 + " AND " +
                        matcher.colName + PatternMatcher.cmp(PatternMatcher.BEFORE, matcher.include2) + "@" + name1 + ")";
                result.params.add(new FuncParam("datetime", name0, matcher.values.get(0)));
                result.params.add(new FuncParam("datetime", name1, matcher.values.get(1)));
                break;
            case PatternMatcher.IS_NULL:
            case PatternMatcher.IS_NOT_NULL:
                result.query = String.format("(%s %s)", matcher.colName, matcher.op);
                break;
            case PatternMatcher.NONE:
            default:
                result.query = "(1 = 1)";
                break;
        }

        return result;
    }

    public String aggrToSql(GroupAggregation aggr) {
        AggrFunctionInfo funcInfo = null;
        for (AggrFunctionInfo info: descriptor.aggregations) {
            if (info.functionName.equals(aggr.aggType)) {
                funcInfo = info;
                break;
            }
        }
        if (funcInfo != null) {
            String brackets = descriptor.nameBrackets;
            String sql = GrokConnectUtil.isNotEmpty(aggr.colName) ? funcInfo.dbFunctionName.replaceAll("#", addBrackets(aggr.colName)) : funcInfo.dbFunctionName;
            return sql + " as " + (GrokConnectUtil.isNotEmpty(aggr.resultColName)
                    ?  brackets.charAt(0) + aggr.resultColName + brackets.charAt(brackets.length() - 1)
                    : brackets.charAt(0) + funcInfo.functionName + "(" + aggr.colName + ")" + brackets.charAt(brackets.length() - 1));
        }
        else
            return null;
    }

    public String addBrackets(String name) {
        String brackets = descriptor.nameBrackets;
        return Arrays.stream(name.split("\\."))
                .map((str) -> str.startsWith(brackets.substring(0, 1)) ? str
                        : brackets.charAt(0) + str + brackets.substring(brackets.length() - 1))
                .collect(Collectors.joining("."));
    }

    public String limitToSql(String query, Integer limit) {
        return query + "limit " + limit.toString();
    }

    public String queryTableSql(TableQuery query) {
        return query.toSql();
    }

    public String insertSql(InsertRows m) {
        if (m.columns == null || m.columns.isEmpty())
            throw new MutationValidationException("InsertRows requires a non-empty columns list");
        m.columns.forEach(this::validateMutationIdentifier);
        return "INSERT INTO " + mutationTableName(m) + " (" +
                m.columns.stream().map(this::addBrackets).collect(Collectors.joining(", ")) +
                ") VALUES (" + String.join(", ", Collections.nCopies(m.columns.size(), "?")) + ")";
    }

    /**
     * INSERT that tolerates duplicate-key conflicts by skipping them (batch mode, {@code errorOnDuplicate ==
     * false}; connector-writes WO-6). The generic provider has no portable ignore-duplicates form, so it
     * emits the plain insert — duplicates surface as errors. Providers with a native form (Postgres
     * {@code ON CONFLICT DO NOTHING}) override; a skipped row reports an update count of 0.
     */
    public String insertIgnoreDuplicatesSql(InsertRows m) {
        return insertSql(m);
    }

    /**
     * Creates the streamed bulk-insert loader for {@code m} on {@code conn} (connector-writes WO-5).
     * The default chunked-{@code executeBatch} loader works on any prepared-statement provider;
     * providers with a native fast path (Postgres COPY) override this.
     */
    public BulkLoader createBulkLoader(Connection conn, InsertRows m) throws SQLException {
        return new BatchInsertBulkLoader(this, conn, m);
    }

    public String updateSql(UpdateRows m, List<FuncParam> collectedParams) {
        if (m.setColumns == null || m.setColumns.isEmpty())
            throw new MutationValidationException("UpdateRows requires a non-empty setColumns list");
        if (m.setValues == null || m.setTypes == null
                || m.setValues.size() != m.setColumns.size() || m.setTypes.size() != m.setColumns.size())
            throw new MutationValidationException("UpdateRows setColumns/setValues/setTypes must be parallel lists of equal size");
        m.setColumns.forEach(this::validateMutationIdentifier);
        StringBuilder sql = new StringBuilder("UPDATE ").append(mutationTableName(m)).append(" SET ");
        sql.append(m.setColumns.stream().map((c) -> addBrackets(c) + " = ?").collect(Collectors.joining(", ")));
        appendWhere(sql, m.whereClauses, m.whereOp, collectedParams);
        return sql.toString();
    }

    /**
     * Key-based batched UPDATE for the {@code mode == "update"} batch path (connector-writes WO-6):
     * {@code UPDATE <fqtn> SET <non-key col> = ? ... WHERE <keyColumn> = ? AND ...}, one row of
     * placeholders bound per CSV / inline row. {@code keyColumns} must be a non-empty subset of
     * {@code columns} (a payload lacking a key column is a validation error) and at least one non-key
     * column must remain to update. Bind parameters follow {@link #updateByKeyBindOrder}.
     */
    public String updateByKeySql(InsertRows m) {
        validateUpdateByKeyColumns(m);
        List<String> nonKey = updateNonKeyColumns(m);
        StringBuilder sql = new StringBuilder("UPDATE ").append(mutationTableName(m)).append(" SET ");
        sql.append(nonKey.stream().map((c) -> addBrackets(c) + " = ?").collect(Collectors.joining(", ")));
        sql.append(" WHERE ").append(m.keyColumns.stream().map((k) -> addBrackets(k) + " = ?")
                .collect(Collectors.joining(" AND ")));
        return sql.toString();
    }

    /** Parameter bind order for {@link #updateByKeySql}: non-key columns (in column order) then key columns. */
    public int[] updateByKeyBindOrder(InsertRows m) {
        validateUpdateByKeyColumns(m);
        List<Integer> order = new ArrayList<>();
        for (int c = 0; c < m.columns.size(); c++)
            if (!m.keyColumns.contains(m.columns.get(c)))
                order.add(c);
        for (String key : m.keyColumns)
            order.add(m.columns.indexOf(key));
        int[] result = new int[order.size()];
        for (int i = 0; i < result.length; i++)
            result[i] = order.get(i);
        return result;
    }

    protected void validateUpdateByKeyColumns(InsertRows m) {
        if (m.columns == null || m.columns.isEmpty())
            throw new MutationValidationException("Update mode requires a non-empty columns list");
        if (m.keyColumns == null || m.keyColumns.isEmpty())
            throw new MutationValidationException("Update mode requires a non-empty keyColumns list");
        m.columns.forEach(this::validateMutationIdentifier);
        for (String key : m.keyColumns) {
            validateMutationIdentifier(key);
            if (!m.columns.contains(key))
                throw new MutationValidationException("keyColumn '" + key + "' must be present in columns");
        }
        if (updateNonKeyColumns(m).isEmpty())
            throw new MutationValidationException("Update mode requires at least one non-key column");
    }

    private List<String> updateNonKeyColumns(InsertRows m) {
        List<String> nonKey = new ArrayList<>();
        for (String col : m.columns)
            if (!m.keyColumns.contains(col))
                nonKey.add(col);
        return nonKey;
    }

    public String deleteSql(DeleteRows m, List<FuncParam> collectedParams) {
        if ((m.whereClauses == null || m.whereClauses.isEmpty()) && !m.allowFullTable)
            throw new MutationValidationException("DeleteRows without whereClauses requires allowFullTable");
        StringBuilder sql = new StringBuilder("DELETE FROM ").append(mutationTableName(m));
        appendWhere(sql, m.whereClauses, m.whereOp, collectedParams);
        return sql.toString();
    }

    /**
     * Dialect-specific INSERT-or-UPDATE for {@code rowCount} value tuples. The generic JDBC provider has
     * no portable upsert — overridden per dialect (connector-writes WO-4); MutationRunner maps the
     * {@link UnsupportedOperationException} to the structured capability error. {@code rowCount} is 1 for
     * dialects executed via addBatch ({@link #upsertBatchRows}==1: Postgres/MySQL/Oracle) and the chunk
     * size for MERGE-over-VALUES dialects (MS SQL, Snowflake).
     */
    public String upsertSql(UpsertRows m, int rowCount) {
        throw new UnsupportedOperationException("Upsert is not supported for provider " + descriptor.type);
    }

    /**
     * Number of value tuples emitted per upsert statement: 1 = single-row statements executed via
     * addBatch (Postgres, MySQL, Oracle); &gt;1 = multi-row VALUES chunks executed one executeUpdate per
     * chunk (MERGE-over-VALUES dialects). {@code columnCount} lets bound-parameter-limited dialects cap
     * the chunk size (MS SQL: 2100-parameter limit).
     */
    public int upsertBatchRows(int columnCount) {
        return 1;
    }

    protected void validateUpsertColumns(UpsertRows m) {
        if (m.columns == null || m.columns.isEmpty())
            throw new MutationValidationException("UpsertRows requires a non-empty columns list");
        if (m.matchKeys == null || m.matchKeys.isEmpty())
            throw new MutationValidationException("UpsertRows requires a non-empty matchKeys list");
        m.columns.forEach(this::validateMutationIdentifier);
        for (String key : m.matchKeys) {
            validateMutationIdentifier(key);
            if (!m.columns.contains(key))
                throw new MutationValidationException("matchKey '" + key + "' must be present in columns");
        }
    }

    protected List<String> upsertNonKeyColumns(UpsertRows m) {
        List<String> nonKey = new ArrayList<>();
        for (String col : m.columns)
            if (!m.matchKeys.contains(col))
                nonKey.add(col);
        return nonKey;
    }

    private String valuesTuples(int columnCount, int rowCount) {
        String tuple = "(" + String.join(", ", Collections.nCopies(columnCount, "?")) + ")";
        return String.join(", ", Collections.nCopies(rowCount, tuple));
    }

    /** MERGE-over-VALUES upsert (MS SQL / Snowflake shape). {@code trailingSemicolon} for T-SQL's MERGE. */
    protected String mergeValuesUpsertSql(UpsertRows m, int rowCount, boolean trailingSemicolon) {
        validateUpsertColumns(m);
        String colList = m.columns.stream().map(this::addBrackets).collect(Collectors.joining(", "));
        String on = m.matchKeys.stream().map((k) -> "t." + addBrackets(k) + " = src." + addBrackets(k))
                .collect(Collectors.joining(" AND "));
        StringBuilder sql = new StringBuilder("MERGE INTO ").append(mutationTableName(m)).append(" AS t USING (VALUES ")
                .append(valuesTuples(m.columns.size(), rowCount)).append(") AS src (").append(colList)
                .append(") ON (").append(on).append(")");
        List<String> nonKey = upsertNonKeyColumns(m);
        if (!nonKey.isEmpty())
            sql.append(" WHEN MATCHED THEN UPDATE SET ").append(nonKey.stream()
                    .map((c) -> "t." + addBrackets(c) + " = src." + addBrackets(c)).collect(Collectors.joining(", ")));
        sql.append(" WHEN NOT MATCHED THEN INSERT (").append(colList).append(") VALUES (")
                .append(m.columns.stream().map((c) -> "src." + addBrackets(c)).collect(Collectors.joining(", ")))
                .append(")");
        if (trailingSemicolon)
            sql.append(";");
        return sql.toString();
    }

    public String mutationTableName(TableMutation m) {
        validateMutationIdentifier(m.tableName);
        if (GrokConnectUtil.isNotEmpty(m.schema))
            validateMutationIdentifier(m.schema);
        if (GrokConnectUtil.isNotEmpty(m.catalog))
            validateMutationIdentifier(m.catalog);
        return SqlNames.fullTableName(m.tableName, m.schema, m.catalog, this);
    }

    /**
     * Neutralizes caller-supplied identifiers on the (destructive) mutation boundary: rejects any
     * table/schema/catalog/column/predicate-field name whose segments are empty or carry the
     * provider's bracket/quote chars or control chars — a hostile identifier becomes a clean
     * pre-validation 400 instead of a downstream db-error. Applied only here, not in the shared
     * addBrackets read path. Dotted names (schema.table, table.column) are validated per segment.
     */
    public void validateMutationIdentifier(String identifier) {
        if (GrokConnectUtil.isEmpty(identifier))
            throw new MutationValidationException("Mutation identifier must not be empty");
        String brackets = descriptor.nameBrackets;
        for (String segment : identifier.split("\\.", -1)) {
            if (segment.isEmpty())
                throw new MutationValidationException("Invalid mutation identifier: '" + identifier + "'");
            for (int i = 0; i < segment.length(); i++) {
                char c = segment.charAt(i);
                if (c < 0x20 || c == '"' || c == '\'' || c == '`' || c == ';' || brackets.indexOf(c) >= 0)
                    throw new MutationValidationException("Illegal character in mutation identifier: '" + identifier + "'");
            }
        }
    }

    /**
     * Offending column name for a per-row mutation error, or {@code null} if the driver does not expose
     * it (connector-writes WO-6). Postgres reads {@code PSQLException.getServerErrorMessage()}; other
     * drivers fall back to code + message only.
     */
    public String mutationErrorColumn(SQLException e) {
        return null;
    }

    /** Human-readable mutation-error message; providers may enrich it (Postgres appends the constraint name). */
    public String mutationErrorMessage(SQLException e) {
        return e.getMessage();
    }

    // ---- DDL emission (connector-writes WO-B6) — statement text only; execution/dryRun is WO-B7 ----

    /** Native column type for a dg scalar type via the descriptor's hand-authored reverse map (ARCHITECTURE §3.3). */
    public String nativeType(String dgType) {
        String nativeType = descriptor.dgToNativeType == null ? null : descriptor.dgToNativeType.get(dgType);
        if (nativeType == null)
            throw new MutationValidationException("No native type for dg type '" + dgType
                    + "' on provider " + descriptor.type);
        return nativeType;
    }

    /**
     * DEFAULT-value literal typed by the dg type (DDL cannot bind parameters): parse-validated numerics,
     * true/false bools, single-quote-doubled strings and ISO datetimes — the DomainDdlGenerator.sqlLiteral
     * discipline. A hostile value either fails validation or stays an inert quoted literal.
     */
    public String ddlLiteral(String dgType, String value) {
        if (value == null)
            return "NULL";
        switch (dgType) {
            case "int":
            case "bigint":
                try {
                    return Long.toString(Long.parseLong(value.trim()));
                } catch (NumberFormatException e) {
                    throw new MutationValidationException("Invalid " + dgType + " literal: '" + value + "'");
                }
            case "float":
                double parsed;
                try {
                    parsed = Double.parseDouble(value.trim());
                } catch (NumberFormatException e) {
                    throw new MutationValidationException("Invalid float literal: '" + value + "'");
                }
                if (Double.isNaN(parsed) || Double.isInfinite(parsed))
                    throw new MutationValidationException("Invalid float literal: '" + value + "'");
                return Double.toString(parsed);
            case "bool":
                return boolDdlLiteral(value.trim().equalsIgnoreCase("true"));
            default: // string, datetime
                return "'" + stringLiteralEscape(value) + "'";
        }
    }

    /**
     * Escapes a value for embedding in a single-quoted SQL literal. Quote doubling is enough for
     * dialects with standard-conforming strings; MySQL/MariaDB treat backslash as an escape character
     * by default and override to double it too.
     */
    protected String stringLiteralEscape(String value) {
        return value.replace("'", "''");
    }

    /** Boolean DDL literal; dialects without a boolean type (MS SQL bit, Oracle number(1)) emit 1/0. */
    protected String boolDdlLiteral(boolean value) {
        return value ? "true" : "false";
    }

    /** Column definition: name, native type, DEFAULT before NOT NULL (the order Oracle requires; valid everywhere). */
    protected String columnDefinitionSql(ColumnSpec c) {
        if (c == null)
            throw new MutationValidationException("Column spec is required");
        validateMutationIdentifier(c.name);
        if (GrokConnectUtil.isEmpty(c.type))
            throw new MutationValidationException("Column '" + c.name + "' requires a dg type");
        StringBuilder def = new StringBuilder(addBrackets(c.name)).append(" ").append(nativeType(c.type));
        if (c.defaultValue != null)
            def.append(" DEFAULT ").append(ddlLiteral(c.type, c.defaultValue));
        if (!c.nullable)
            def.append(" NOT NULL");
        return def.toString();
    }

    /** IF NOT EXISTS on CREATE TABLE; dialects without the form (MS SQL, Oracle) emit a plain CREATE — a duplicate table surfaces as a db-error. */
    protected boolean supportsCreateIfNotExists() {
        return true;
    }

    /** CREATE TABLE (+ one CREATE [UNIQUE] INDEX per {@link IndexSpec}); FKs are added afterwards via {@link AddKey}. */
    public List<String> createTableSql(CreateTable m) {
        if (m.columns == null || m.columns.isEmpty())
            throw new MutationValidationException("CreateTable requires a non-empty columns list");
        String table = mutationTableName(m);
        List<String> defs = new ArrayList<>();
        for (ColumnSpec c : m.columns)
            defs.add(columnDefinitionSql(c));
        if (m.primaryKey != null && !m.primaryKey.isEmpty()) {
            m.primaryKey.forEach(this::validateMutationIdentifier);
            defs.add("PRIMARY KEY (" + m.primaryKey.stream().map(this::addBrackets)
                    .collect(Collectors.joining(", ")) + ")");
        }
        List<String> statements = new ArrayList<>();
        statements.add("CREATE TABLE " + (m.ifNotExists && supportsCreateIfNotExists() ? "IF NOT EXISTS " : "")
                + table + " (" + String.join(", ", defs) + ")");
        if (m.indexes != null)
            for (IndexSpec index : m.indexes)
                statements.add(createIndexStatement(m, index.name, index.columns, index.unique));
        return statements;
    }

    /** Single-action ALTER TABLE (ARCHITECTURE §3.4.1); per-action field presence is validated here. */
    public String alterTableSql(AlterTable m) {
        String table = mutationTableName(m);
        if (GrokConnectUtil.isEmpty(m.action))
            throw new MutationValidationException("AlterTable requires an action");
        switch (m.action) {
            case "addColumn":
                if (m.column == null)
                    throw new MutationValidationException("AlterTable addColumn requires a column spec");
                return alterAddColumnSql(m, table);
            case "dropColumn":
                validateMutationIdentifier(m.columnName);
                return "ALTER TABLE " + table + " DROP COLUMN " + addBrackets(m.columnName);
            case "renameColumn":
                validateMutationIdentifier(m.columnName);
                validateMutationIdentifier(m.newName);
                return alterRenameColumnSql(m, table);
            case "changeType":
                validateMutationIdentifier(m.columnName);
                if (GrokConnectUtil.isEmpty(m.newType))
                    throw new MutationValidationException("AlterTable changeType requires newType");
                return alterChangeTypeSql(m, table);
            case "setNullable":
                validateMutationIdentifier(m.columnName);
                if (m.nullable == null)
                    throw new MutationValidationException("AlterTable setNullable requires an explicit nullable value");
                return alterSetNullableSql(m, table);
            default:
                throw new MutationValidationException("Unknown AlterTable action: '" + m.action + "'");
        }
    }

    protected String alterAddColumnSql(AlterTable m, String table) {
        return "ALTER TABLE " + table + " ADD COLUMN " + columnDefinitionSql(m.column);
    }

    protected String alterRenameColumnSql(AlterTable m, String table) {
        return "ALTER TABLE " + table + " RENAME COLUMN " + addBrackets(m.columnName)
                + " TO " + addBrackets(m.newName);
    }

    /** PG shape; MySQL (MODIFY), MS SQL and Oracle override. */
    protected String alterChangeTypeSql(AlterTable m, String table) {
        return "ALTER TABLE " + table + " ALTER COLUMN " + addBrackets(m.columnName)
                + " TYPE " + nativeType(m.newType);
    }

    /** PG shape; dialects that restate the column type (MySQL, MS SQL) or use MODIFY (Oracle) override. */
    protected String alterSetNullableSql(AlterTable m, String table) {
        return "ALTER TABLE " + table + " ALTER COLUMN " + addBrackets(m.columnName)
                + (m.nullable ? " DROP NOT NULL" : " SET NOT NULL");
    }

    public String createIndexSql(CreateIndex m) {
        return createIndexStatement(m, m.indexName, m.columns, m.unique);
    }

    /** CREATE [UNIQUE] INDEX; auto-name ix_&lt;table&gt;_&lt;col1&gt;_.. (the domain generator convention) when absent. */
    protected String createIndexStatement(TableMutation m, String name, List<String> columns, boolean unique) {
        if (columns == null || columns.isEmpty())
            throw new MutationValidationException("Index requires a non-empty columns list");
        columns.forEach(this::validateMutationIdentifier);
        String indexName = GrokConnectUtil.isNotEmpty(name) ? name
                : "ix_" + baseTableName(m) + "_" + String.join("_", columns);
        validateMutationIdentifier(indexName);
        return "CREATE " + (unique ? "UNIQUE " : "") + "INDEX " + addBrackets(indexName)
                + " ON " + mutationTableName(m) + " ("
                + columns.stream().map(this::addBrackets).collect(Collectors.joining(", ")) + ")";
    }

    /** Unqualified table name for auto-generated index/constraint names. */
    protected String baseTableName(TableMutation m) {
        validateMutationIdentifier(m.tableName);
        return m.tableName.substring(m.tableName.lastIndexOf('.') + 1);
    }

    /** PG/Oracle shape (schema-qualified name); MySQL / MS SQL override with DROP INDEX &lt;name&gt; ON &lt;table&gt;. */
    public String dropIndexSql(DropIndex m) {
        validateMutationIdentifier(m.indexName);
        String qualified = m.indexName;
        if (GrokConnectUtil.isNotEmpty(m.schema)) {
            validateMutationIdentifier(m.schema);
            qualified = m.schema + "." + m.indexName;
        }
        return "DROP INDEX " + addBrackets(qualified);
    }

    /** ALTER TABLE ADD CONSTRAINT — primary or foreign key; auto-names pk_&lt;table&gt; / fk_&lt;table&gt;_&lt;col&gt; (the domain convention). */
    public String addKeySql(AddKey m) {
        String table = mutationTableName(m);
        if (m.columns == null || m.columns.isEmpty())
            throw new MutationValidationException("AddKey requires a non-empty columns list");
        m.columns.forEach(this::validateMutationIdentifier);
        String colList = m.columns.stream().map(this::addBrackets).collect(Collectors.joining(", "));
        if ("primary".equals(m.keyType)) {
            String name = GrokConnectUtil.isNotEmpty(m.constraintName) ? m.constraintName : "pk_" + baseTableName(m);
            validateMutationIdentifier(name);
            return "ALTER TABLE " + table + " ADD CONSTRAINT " + addBrackets(name) + " PRIMARY KEY (" + colList + ")";
        }
        if ("foreign".equals(m.keyType)) {
            if (GrokConnectUtil.isEmpty(m.refTable) || m.refColumns == null || m.refColumns.isEmpty())
                throw new MutationValidationException("AddKey foreign requires refTable and refColumns");
            validateMutationIdentifier(m.refTable);
            m.refColumns.forEach(this::validateMutationIdentifier);
            String name = GrokConnectUtil.isNotEmpty(m.constraintName) ? m.constraintName
                    : "fk_" + baseTableName(m) + "_" + String.join("_", m.columns);
            validateMutationIdentifier(name);
            return "ALTER TABLE " + table + " ADD CONSTRAINT " + addBrackets(name) + " FOREIGN KEY (" + colList
                    + ") REFERENCES " + addBrackets(m.refTable) + " ("
                    + m.refColumns.stream().map(this::addBrackets).collect(Collectors.joining(", ")) + ")";
        }
        throw new MutationValidationException("AddKey keyType must be 'primary' or 'foreign', got '" + m.keyType + "'");
    }

    public String dropTableSql(DropTable m) {
        return "DROP TABLE " + mutationTableName(m);
    }

    public String truncateTableSql(TruncateTable m) {
        return "TRUNCATE TABLE " + mutationTableName(m);
    }

    /** Appends a WHERE section identical in shape to TableQuery.toSql's (shared PredicateCompiler). */
    protected void appendWhere(StringBuilder sql, List<FieldPredicate> whereClauses, String whereOp, List<FuncParam> collectedParams) {
        if (whereClauses == null || whereClauses.isEmpty())
            return;
        String op = GrokConnectUtil.isEmpty(whereOp) ? "and" : whereOp;
        if (!op.equals("and") && !op.equals("or"))
            throw new MutationValidationException("whereOp must be 'and' or 'or', got '" + whereOp + "'");
        List<String> clauses = new ArrayList<>();
        for (FieldPredicate clause : whereClauses) {
            validateMutationIdentifier(clause.field);
            clauses.add(String.format("  (%s)", PredicateCompiler.compile(clause, this, collectedParams)));
        }
        sql.append(System.lineSeparator()).append("WHERE").append(System.lineSeparator());
        sql.append(String.join(String.format(" %s%s", op, System.lineSeparator()), clauses));
    }

    public String castParamValueToSqlDateTime(FuncParam param) {
        return "datetime('" + param.value.toString() + "')";
    }

    public java.util.Properties defaultConnectionProperties(DataConnection conn) {
        java.util.Properties properties = new java.util.Properties();
        if (conn.credentials != null) {
            String method = (String) conn.credentials.parameters.get("#chosen-auth-method");
            boolean includeLogin = method == null || descriptor == null
                    || descriptor.credentialsTemplate == null
                    || descriptor.credentialsTemplate.stream()
                        .filter(p -> p.name.equals(DbCredentials.LOGIN))
                        .anyMatch(p -> p.category != null && p.category.contains(method));
            if (includeLogin) {
                setIfNotNull(properties, "user", conn.credentials.getLogin());
                setIfNotNull(properties, "password", conn.credentials.getPassword());
            }
        }
        return properties;
    }

    public ResultSetManager getResultSetManager() {
        return DefaultResultSetManager.getDefaultManager();
    }

    public static void setIfNotNull(java.util.Properties properties, String key, String value) {
        if (value != null)
            properties.setProperty(key, value);
    }

    public static void setIfNotEmpty(java.util.Properties properties, String key, String value) {
        if (GrokConnectUtil.isNotEmpty(value))
            properties.setProperty(key, value);
    }
}

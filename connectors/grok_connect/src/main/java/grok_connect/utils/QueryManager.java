package grok_connect.utils;

import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import grok_connect.GrokConnect;
import grok_connect.connectors_info.FuncCall;
import grok_connect.log.QueryLogger;
import grok_connect.log.QueryLoggerImpl;
import grok_connect.providers.JdbcDataProvider;
import org.slf4j.Logger;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import serialization.Column;
import serialization.DataFrame;

public class QueryManager {
    private static final Gson gson = new GsonBuilder()
            .registerTypeAdapter(Property.class, new PropertyAdapter())
            .create();
    private static final String FETCH_SIZE_KEY = "connectFetchSize";
    private static final int MAX_CHUNK_SIZE_BYTES = 10_000_000;
    private static final int MAX_FETCH_SIZE = 100000;
    private static final int MIN_FETCH_SIZE = 100;
    private static final int FIRST_FETCH_SIZE = 100;
    private final JdbcDataProvider provider;
    private final FuncCall query;
    private final QueryLogger queryLogger;
    private final Logger logger;
    private int currentFetchSize = FIRST_FETCH_SIZE;
    private int chunkSize = -1;
    private SchemeInfo schemeInfo;
    private ResultSet resultSet;
    private Connection connection;
    private boolean changedFetchSize;
    private boolean supportTransactions;

    public QueryManager(String message) {
        query = gson.fromJson(message, FuncCall.class);
        queryLogger = new QueryLoggerImpl(query.printLevels);
        logger = queryLogger.getLogger();
        query.log = "";
        query.setParamValues();
        query.afterDeserialization();
        provider = GrokConnect.providerManager.getByName(query.func.connection.dataSource);
        if (query.func.options != null && query.func.options.containsKey(FETCH_SIZE_KEY)) {
            String value = query.func.options.get(FETCH_SIZE_KEY).toString();
            Pattern pattern = Pattern.compile("(\\d+) MB");
            Matcher matcher = pattern.matcher(value);
            if (matcher.find()) {
                chunkSize = Integer.parseInt(matcher.group(1)) * 1_000_000;
            } else {
                changedFetchSize = true;
                currentFetchSize = Integer.parseInt(value);
            }
        }
    }

    public void initResultSet() throws ClassNotFoundException, GrokConnectException, QueryCancelledByUser, SQLException {
        logger.debug("Initializing resultSet");
        long startTime = System.currentTimeMillis();
        connection = provider.getConnection(query.func.connection);
        resultSet = provider.getResultSet(query, connection, logger);
        supportTransactions = connection.getMetaData().supportsTransactions();
        logger.debug("Finished resultSet init, execution time: {} ms", System.currentTimeMillis() - startTime);
    }

    public void initScheme() throws QueryCancelledByUser, SQLException {
        logger.debug("Initializing schemeInfo");
        long startTime = System.currentTimeMillis();
        if (resultSet == null) {
            logger.debug("resultSet is null, return empty schemeInfo");
            schemeInfo = new SchemeInfo(new ArrayList<>(), new ArrayList<>(), new ArrayList<>());
            return;
        }
        schemeInfo = provider.resultSetScheme(query, resultSet, logger);
        logger.debug("Finished schemeInfo init, execution time: {} ms", System.currentTimeMillis() - startTime);
    }

    public DataFrame getSubDF() throws IOException, SQLException, QueryCancelledByUser {
        logger.trace("getSubDF was called with argument");
        List<Column> columns = schemeInfo.columns;
        for (Column column : columns) {
            column.empty();
        }
        DataFrame df = new DataFrame();
        if (!connection.isClosed() && !resultSet.isClosed()) {
            logger.debug("Receiving part of resultSet");
            long startTime = System.currentTimeMillis();
            df = provider.getResultSetSubDf(query, resultSet, columns,
            schemeInfo.supportedType, schemeInfo.initColumn, currentFetchSize, logger);
            logger.debug("Received part of resultSet, execution time: {} ms", System.currentTimeMillis() - startTime);
        }

        if (supportTransactions && df.rowCount != 0 && !changedFetchSize) {
            currentFetchSize = getFetchSize(df);
            logger.debug("Fetch size: {}", currentFetchSize);
            if (!provider.descriptor.type.equals("Virtuoso")) {
                resultSet.setFetchSize(currentFetchSize);
            }
        }
        return df;
    }

    public void closeConnection() throws SQLException {
        if (connection != null && !connection.isClosed()) {
            if (!connection.getAutoCommit())
                connection.commit();
            provider.providerManager.getQueryMonitor().removeResultSet(query.id);
            connection.close();
        } else {
            provider.providerManager.getQueryMonitor().removeResultSet(query.id);
        }
    }

    public boolean isResultSetInitialized() {
        return resultSet != null;
    }

    public FuncCall getQuery() {
        return query;
    }

    public QueryLogger getQueryLogger() {
        return queryLogger;
    }

    private int getFetchSize(DataFrame df) {
        int maxChunkSize = chunkSize != -1 ? chunkSize : MAX_CHUNK_SIZE_BYTES;
        int fetchSize = Math.round(maxChunkSize / (float) (df.memoryInBytes() / df.rowCount));
        return Math.min(Math.max(MIN_FETCH_SIZE, fetchSize), MAX_FETCH_SIZE);
    }
}

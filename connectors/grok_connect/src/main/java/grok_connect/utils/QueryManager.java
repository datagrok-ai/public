package grok_connect.utils;

import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import grok_connect.GrokConnect;
import grok_connect.connectors_info.FuncCall;
import grok_connect.log.QueryLogger;
import grok_connect.log.QueryLoggerImpl;
import grok_connect.providers.JdbcDataProvider;
import org.slf4j.Logger;
import serialization.*;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class QueryManager {
    private static final Gson gson = new GsonBuilder()
            .registerTypeAdapter(Property.class, new PropertyAdapter())
            .create();
    private static final String AUX_FETCH_SIZE_KEY = "fetchSize";
    private static final int DEFAULT_FETCH_SIZE = 10000;
    private final JdbcDataProvider provider;
    private final FuncCall query;
    private final QueryLogger queryLogger;
    private final Logger logger;
    private SchemeInfo schemeInfo;
    private ResultSet resultSet;
    private Connection connection;
    private boolean changedFetchSize;

    public QueryManager(String message) {
        query = gson.fromJson(message, FuncCall.class);
        queryLogger = new QueryLoggerImpl(query.printLevels);
        logger = queryLogger.getLogger();
        query.log = "";
        query.setParamValues();
        query.afterDeserialization();
        provider = GrokConnect.providerManager.getByName(query.func.connection.dataSource);
    }

    public void initResultSet() throws ClassNotFoundException, GrokConnectException, QueryCancelledByUser, SQLException {
        logger.debug("Initializing resultSet");
        long startTime = System.currentTimeMillis();
        connection = provider.getConnection(query.func.connection);
        resultSet = provider.getResultSet(query, connection, logger);
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

    public DataFrame getSubDF(int maxIterations) throws IOException, SQLException, QueryCancelledByUser {
        logger.trace("getSubDF was called with argument: {}", maxIterations);
        List<Column> columns = schemeInfo.columns;
        for (Column column : columns) {
            column.empty();
        }
        DataFrame df = new DataFrame();
        if (!connection.isClosed() && !resultSet.isClosed()) {
            logger.debug("Receiving part of resultSet");
            long startTime = System.currentTimeMillis();
            df = provider.getResultSetSubDf(query, resultSet, columns,
            schemeInfo.supportedType, schemeInfo.initColumn, maxIterations, logger);
            logger.debug("Received part of resultSet, execution time: {} ms", System.currentTimeMillis() - startTime);
        }

        if (connection.getMetaData().supportsTransactions() && df.rowCount != 0) {
            int dataFrameWeight = df.memoryInBytes();
            Double memInBytes = dataFrameWeight / 1024.0;
            double fetchSize = df.rowCount / memInBytes * 30000;
            logger.debug("Dataframe weight: {} bytes", dataFrameWeight);
            if (fetchSize > 40000)
                fetchSize = 40000;
            logger.debug("Fetch size: {}", fetchSize);
            if (!provider.descriptor.type.equals("Virtuoso")) {
                resultSet.setFetchSize((int)fetchSize);
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
}

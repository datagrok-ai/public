package grok_connect.utils;

import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.time.LocalDateTime;
import java.time.ZoneOffset;
import java.util.ArrayList;
import java.util.List;
import grok_connect.GrokConnect;
import grok_connect.connectors_info.FuncCall;
import grok_connect.handlers.SessionHandler;
import grok_connect.providers.JdbcDataProvider;
import org.apache.hadoop.yarn.webapp.hamlet2.Hamlet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.*;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class QueryManager {
    private static final Logger LOGGER = LoggerFactory.getLogger(QueryManager.class);
    private static final Gson gson = new GsonBuilder()
            .registerTypeAdapter(Property.class, new PropertyAdapter())
            .create();
    private static final String AUX_FETCH_SIZE_KEY = "fetchSize";
    private static final int DEFAULT_FETCH_SIZE = 10000;
    private final JdbcDataProvider provider;
    private final FuncCall query;
    private SchemeInfo schemeInfo;
    private ResultSet resultSet;
    private Connection connection;
    private boolean changedFetchSize;

    public QueryManager(FuncCall query, JdbcDataProvider provider) {
        this.query = query;
        this.provider = provider;
    }

    public QueryManager(String message) {
        query = gson.fromJson(message, FuncCall.class);
        query.log = "";
        query.setParamValues();
        query.afterDeserialization();
        LOGGER.debug("Initializing with query: {}", query.func.query);
        // DateTime startTime = DateTime.now();
        provider = GrokConnect.providerManager.getByName(query.func.connection.dataSource);
    }

    public void initResultSet() throws ClassNotFoundException, GrokConnectException, QueryCancelledByUser, SQLException {
        LOGGER.trace("initResultSet was called");
        LocalDateTime startTime = LocalDateTime.now();
        if (query.debugQuery) {
            query.log += String.format(SessionHandler.LOG_MESSAGE, String.format("Initializing resultSet, current time: %s",
                    startTime));
        }
        connection = provider.getConnection(query.func.connection);
        resultSet = provider.getResultSet(query, connection);
        if (query.debugQuery) {
            query.log += String.format(SessionHandler.LOG_MESSAGE, String.format("Finished resultSet init, execution time: %s",
                    LocalDateTime.now().toInstant(ZoneOffset.UTC).toEpochMilli()
                            - startTime.toInstant(ZoneOffset.UTC).toEpochMilli()));
        }
        LOGGER.trace("ResultSet received");
    }

    public void initScheme() throws QueryCancelledByUser, SQLException {
        LocalDateTime startTime = LocalDateTime.now();
        if (query.debugQuery) {
            query.log += String.format(SessionHandler.LOG_MESSAGE, String.format("Initializing schemeInfo, current time: %s",
                    startTime));
        }
        if (resultSet == null) {
            schemeInfo = new SchemeInfo(new ArrayList<>(), new ArrayList<>(), new ArrayList<>());
            return;
        }
        schemeInfo = provider.resultSetScheme(query, resultSet);
        if (query.debugQuery) {
            query.log += String.format(SessionHandler.LOG_MESSAGE, String.format("Finished schemeInfo init, execution time: %s",
                    LocalDateTime.now().toInstant(ZoneOffset.UTC).toEpochMilli()
                            - startTime.toInstant(ZoneOffset.UTC).toEpochMilli()));
        }
        LOGGER.trace("SchemeInfo initialization finished");
    }

    public DataFrame getSubDF(int maxIterations) throws IOException, SQLException, QueryCancelledByUser {
        LOGGER.trace("getSubDF was called with argument: {}", maxIterations);
        List<Column> columns = schemeInfo.columns;
        for (Column column : columns) {
            column.empty();
        }
        DataFrame df = new DataFrame();
        if (!connection.isClosed() && !resultSet.isClosed()) {
            LocalDateTime startTime = LocalDateTime.now();
            LOGGER.trace("Calling getResultSetSubDf");
            if (query.debugQuery) {
                query.log += String.format(SessionHandler.LOG_MESSAGE, String.format("Receiving part of resultSet, "
                                + "current time: %s",
                        startTime));
            }
            df = provider.getResultSetSubDf(query, resultSet, columns,
            schemeInfo.supportedType, schemeInfo.initColumn, maxIterations);
            if (query.debugQuery) {
                query.log += String.format(SessionHandler.LOG_MESSAGE, String.format("Received part of resultSet, execution time: %s",
                        LocalDateTime.now().toInstant(ZoneOffset.UTC).toEpochMilli()
                                - startTime.toInstant(ZoneOffset.UTC).toEpochMilli()));
            }
        }

        if (connection.getMetaData().supportsTransactions() && df.rowCount != 0) {
            Double memInBytes = df.memoryInBytes() / 1024.0;
            double fetchSize = df.rowCount / memInBytes * 30000;
            if (query.debugQuery) {
                query.log += String.format(SessionHandler.LOG_MESSAGE, String.format("Dataframe weight %s",
                        memInBytes));
            }
            if (fetchSize > 40000)
                fetchSize = 40000;
            query.log += String.format(SessionHandler.LOG_MESSAGE, String.format("New fetch size: %s", fetchSize));
            if (!provider.descriptor.type.equals("Virtuoso")) {
                resultSet.setFetchSize((int)fetchSize);
                LOGGER.trace("Set new fetch size to {}", fetchSize);
            }
        }
        return df;
    }

    public void closeConnection() throws SQLException {
        LOGGER.trace("Closing connection");
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
}

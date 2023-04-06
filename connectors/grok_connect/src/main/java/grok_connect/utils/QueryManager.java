package grok_connect.utils;

import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import grok_connect.GrokConnect;
import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.JdbcDataProvider;
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
        connection = provider.getConnection(query.func.connection);
        resultSet = provider.getResultSet(query, connection);
        LOGGER.trace("ResultSet received");
    }

    public void initScheme() throws QueryCancelledByUser, SQLException {
        LOGGER.trace("Starting schemeInfo initialization");
        if (resultSet == null) {
            schemeInfo = new SchemeInfo(new ArrayList<>(), new ArrayList<>(), new ArrayList<>());
            return;
        }
        schemeInfo = provider.resultSetScheme(query, resultSet);
        LOGGER.trace("SchemeInfo initialization finished");
    }

    public DataFrame getSubDF(int maxIterations) throws IOException, SQLException, QueryCancelledByUser {
        LOGGER.trace("getSubDF was called with argument: {}", maxIterations);
        List<Column> columns = schemeInfo.columns;
        for (Column column : columns) {
            column.empty();
        }
        if (!changedFetchSize && (!query.aux.containsKey(AUX_FETCH_SIZE_KEY)
                || query.aux.get(AUX_FETCH_SIZE_KEY).equals("dynamic"))) {
            if (connection.getMetaData().supportsTransactions())
                resultSet.setFetchSize(DEFAULT_FETCH_SIZE);
            LOGGER.trace("Changing fetchSize");
            changedFetchSize = true;
        }
        DataFrame df = new DataFrame();
        if (!connection.isClosed() && !resultSet.isClosed()) {
            LOGGER.trace("Calling getResultSetSubDf");
            df = provider.getResultSetSubDf(query, resultSet, columns,
                schemeInfo.supportedType, schemeInfo.initColumn, maxIterations);
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

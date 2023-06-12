package grok_connect.utils;

import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import grok_connect.GrokConnect;
import grok_connect.connectors_info.FuncCall;
import grok_connect.log.EventType;
import grok_connect.providers.JdbcDataProvider;
import org.slf4j.Logger;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import org.slf4j.Marker;
import serialization.Column;
import serialization.DataFrame;

public class QueryManager {
    public static final String CHUNK_NUMBER_TAG = "chunkNumber";
    private static final Gson gson = new GsonBuilder()
            .registerTypeAdapter(Property.class, new PropertyAdapter())
            .create();
    public static final String FETCH_SIZE_KEY = "connectFetchSize";
    public static final String DRY_RUN_KEY = "dryRun";
    private static final int MAX_CHUNK_SIZE_BYTES = 10_000_000;
    private static final int MAX_FETCH_SIZE = 100000;
    private static final int MIN_FETCH_SIZE = 100;
    public boolean isDryRun = false;
    private final JdbcDataProvider provider;
    private final FuncCall query;
    private final Logger logger;
    private int currentFetchSize = MIN_FETCH_SIZE;
    private int chunkSize = -1;
    private SchemeInfo schemeInfo;
    private ResultSet resultSet;
    private Connection connection;
    private boolean changedFetchSize;
    private boolean supportTransactions;

    public QueryManager(String message, Logger logger) {
        query = gson.fromJson(message, FuncCall.class);
        query.setParamValues();
        query.afterDeserialization();
        this.logger = logger;
        provider = GrokConnect.providerManager.getByName(query.func.connection.dataSource);
        if (query.func.options != null) {
            if (query.func.options.containsKey(FETCH_SIZE_KEY)) {
                setInitFetchSize(query.func.options.get(FETCH_SIZE_KEY).toString());
            }
            if (query.func.options.containsKey(DRY_RUN_KEY)) {
                isDryRun = Boolean.parseBoolean(query.func.options.get(DRY_RUN_KEY).toString());
            }
        }
        if (query.func.aux != null) {
            if (query.func.aux.containsKey(FETCH_SIZE_KEY)) {
                setInitFetchSize(query.func.aux.get(FETCH_SIZE_KEY).toString());
            }
            if (query.func.aux.containsKey(DRY_RUN_KEY)) {
                isDryRun = (Boolean) query.func.aux.get(DRY_RUN_KEY);
            }
        }
    }

    public void initResultSet() throws ClassNotFoundException, GrokConnectException, QueryCancelledByUser, SQLException {
        logger.debug(EventType.RESULT_SET_INIT.getMarker(EventType.Stage.START), "Initializing resultSet");
        logger.debug(EventType.CONNECTION_RECEIVING.getMarker(EventType.Stage.START), "Receiving connection to db");
        connection = provider.getConnection(query.func.connection);
        logger.debug(EventType.CONNECTION_RECEIVING.getMarker(EventType.Stage.END), "Connection was received");
        resultSet = provider.getResultSet(query, connection, logger, currentFetchSize);
        logger.debug(EventType.RESULT_SET_INIT.getMarker(EventType.Stage.END), "Finished resultSet init");
        supportTransactions = connection.getMetaData().supportsTransactions();
    }

    public void initScheme() throws QueryCancelledByUser, SQLException {
        logger.debug(EventType.SCHEME_INFO_INIT.getMarker(EventType.Stage.START), "Initializing schemeInfo");
        if (resultSet == null) {
            logger.debug(EventType.SCHEME_INFO_INIT.getMarker(EventType.Stage.END), "resultSet is null, return empty schemeInfo");
            schemeInfo = new SchemeInfo(new ArrayList<>(), new ArrayList<>(), new ArrayList<>());
            return;
        }
        schemeInfo = provider.resultSetScheme(query, resultSet, logger);
        logger.debug(EventType.SCHEME_INFO_INIT.getMarker(EventType.Stage.END), "Finished schemeInfo init");
    }

    public DataFrame getSubDF(int dfNumber) throws IOException, SQLException, QueryCancelledByUser {
        logger.trace(EventType.MISC.getMarker(), "getSubDF was called with argument");
        logger.debug(EventType.DATAFRAME_PROCESSING.getMarker(dfNumber, EventType.Stage.START), "DataFrame processing was started");
        List<Column> columns = schemeInfo.columns;
        for (Column column : columns) {
            column.empty();
        }
        DataFrame df = new DataFrame();

        if (!connection.isClosed() && !resultSet.isClosed()) {
            df = getResultSetSubDf(dfNumber, columns);
            if (isDryRun) {
                serializeDf(df, dfNumber);
                logger.debug(EventType.DATAFRAME_PROCESSING.getMarker(dfNumber, EventType.Stage.END), "DataFrame processing was finished");
            }
        }

        changeFetchSize(df, dfNumber);

        if (isDryRun && df.rowCount != 0) {
            while (df.rowCount != 0) {
                logger.debug(EventType.DATAFRAME_PROCESSING.getMarker(++dfNumber, EventType.Stage.START), "DataFrame processing was started");
                for (Column column : columns) {
                    column.empty();
                }
                df = getResultSetSubDf(dfNumber, columns);
                serializeDf(df, dfNumber);
                logger.debug(EventType.DATAFRAME_PROCESSING.getMarker(dfNumber, EventType.Stage.END), "DataFrame processing was finished");
                changeFetchSize(df, dfNumber);
            }
        }

        df.tags = new LinkedHashMap<>();
        df.tags.put(CHUNK_NUMBER_TAG, String.valueOf(dfNumber));
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

    private DataFrame getResultSetSubDf(int dfNumber, List<Column> columns) throws QueryCancelledByUser, SQLException, IOException {
        Marker start = EventType.RESULT_SET_PROCESSING.getMarker(dfNumber, EventType.Stage.START);
        Marker end = EventType.RESULT_SET_PROCESSING.getMarker(dfNumber, EventType.Stage.END);
        logger.debug(start, "ResultSet processing was started");
        DataFrame df = provider.getResultSetSubDf(query, resultSet, columns,
                schemeInfo.supportedType, schemeInfo.initColumn, currentFetchSize, logger, dfNumber);
        logger.debug(end, "ResultSet processing has finished");
        return df;
    }

    private void changeFetchSize(DataFrame df, int dfNumber) throws SQLException {
        if (supportTransactions && df.rowCount != 0 && !changedFetchSize) {
            currentFetchSize = getFetchSize(df);
            logger.debug(EventType.MISC.getMarker(dfNumber), "Fetch size: {}", currentFetchSize);
            if (!provider.descriptor.type.equals("Virtuoso")) {
                resultSet.setFetchSize(currentFetchSize);
            }
        }
    }

    private int getFetchSize(DataFrame df) {
        int maxChunkSize = chunkSize != -1 ? chunkSize : MAX_CHUNK_SIZE_BYTES;
        int fetchSize = Math.round(maxChunkSize / (float) (df.memoryInBytes() / df.rowCount));
        return Math.min(Math.max(MIN_FETCH_SIZE, fetchSize), MAX_FETCH_SIZE);
    }

    private void setInitFetchSize(String optionValue) {
        logger.debug(EventType.MISC.getMarker(), "Setting manual fetch size {}", optionValue);
        Pattern pattern = Pattern.compile("(\\d+) MB");
        Matcher matcher = pattern.matcher(optionValue);
        if (matcher.find()) {
            chunkSize = Integer.parseInt(matcher.group(1)) * 1_000_000;
        } else {
            changedFetchSize = true;
            currentFetchSize = Integer.parseInt(optionValue);
        }
    }

    private byte[] serializeDf(DataFrame df, int dfNumber) {
        Marker start = EventType.DATAFRAME_TO_BYTEARRAY_CONVERTING.getMarker(dfNumber, EventType.Stage.START);
        Marker finish = EventType.DATAFRAME_TO_BYTEARRAY_CONVERTING.getMarker(dfNumber, EventType.Stage.END);
        logger.debug(start, "Converting dataframe to byteArray");
        byte[] bytes = df.toByteArray();
        logger.debug(finish, "DataFrame with id {} was converted to byteArray", dfNumber);
        return bytes;
    }
}

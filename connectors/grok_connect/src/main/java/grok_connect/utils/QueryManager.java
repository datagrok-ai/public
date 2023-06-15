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
import serialization.Column;
import serialization.DataFrame;

public class QueryManager {
    public static final String CHUNK_NUMBER_TAG = "chunkNumber";
    private static final Gson gson = new GsonBuilder()
            .registerTypeAdapter(Property.class, new PropertyAdapter())
            .create();
    public static final String FETCH_SIZE_KEY = "connectFetchSize";
    public static final String INIT_FETCH_SIZE_KEY = "initConnectFetchSize";
    public static final String DRY_RUN_KEY = "dryRun";
    private static final int MAX_CHUNK_SIZE_BYTES = 10_000_000;
    private static final int MAX_FETCH_SIZE = 100000;
    private static final int MIN_FETCH_SIZE = 100;
    public boolean dryRun = false;
    private final JdbcDataProvider provider;
    private final FuncCall query;
    private final Logger logger;
    private int currentFetchSize = MIN_FETCH_SIZE;
    private int initFetchSize = MIN_FETCH_SIZE;
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
                setFetchSize(query.func.options.get(FETCH_SIZE_KEY).toString());
            }
            if (query.func.options.containsKey(DRY_RUN_KEY)) {
                dryRun = Boolean.parseBoolean(query.func.options.get(DRY_RUN_KEY).toString());
            }
            if (query.func.options.containsKey(INIT_FETCH_SIZE_KEY)) {
                setInitFetchSize(query.func.options.get(INIT_FETCH_SIZE_KEY).toString());
            }
        }
        if (query.func.aux != null) {
            if (query.func.aux.containsKey(FETCH_SIZE_KEY)) {
                setFetchSize(query.func.aux.get(FETCH_SIZE_KEY).toString());
            }
            if (query.func.aux.containsKey(DRY_RUN_KEY)) {
                dryRun = (Boolean) query.func.aux.get(DRY_RUN_KEY);
            }
            if (query.func.aux.containsKey(INIT_FETCH_SIZE_KEY)) {
                setInitFetchSize(query.func.aux.get(INIT_FETCH_SIZE_KEY).toString());
            }
        }
    }

    public void initResultSet() throws ClassNotFoundException, GrokConnectException, QueryCancelledByUser, SQLException {
        logger.debug(EventType.CONNECTION_RECEIVING.getMarker(EventType.Stage.START), "Receiving connection to db");
        connection = provider.getConnection(query.func.connection);
        logger.debug(EventType.CONNECTION_RECEIVING.getMarker(EventType.Stage.END), "Connection was received");
        resultSet = provider.getResultSet(query, connection, logger, initFetchSize);
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
        DataFrame df = new DataFrame();
        if (resultSet.isAfterLast() && !resultSet.isBeforeFirst()) {
            return df;
        }
        logger.trace(EventType.MISC.getMarker(), "getSubDF was called with argument");
        List<Column> columns = schemeInfo.columns;
        for (Column column : columns) {
            column.empty();
        }

        if (!connection.isClosed() && !resultSet.isClosed()) {
            df = getResultSetSubDf(dfNumber, columns);
        }

        changeFetchSize(df, dfNumber);

        if (dryRun) {
            while (!resultSet.isAfterLast()) {
                for (Column column : columns) {
                    column.empty();
                }
                df = getResultSetSubDf(++dfNumber, columns);
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
        int rowsNumber = dfNumber == 1 ? initFetchSize : currentFetchSize;
        return provider.getResultSetSubDf(query, resultSet, columns,
                schemeInfo.supportedType, schemeInfo.initColumn, rowsNumber, logger, dfNumber, dryRun);
    }

    private void changeFetchSize(DataFrame df, int dfNumber) throws SQLException {
        if (supportTransactions && dryRun || supportTransactions && changedFetchSize) {
            resultSet.setFetchSize(currentFetchSize);
            return;
        }
        if (supportTransactions && df.rowCount != 0) {
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

    private void setFetchSize(String optionValue) {
        if (optionValue.isEmpty()) {
            return;
        }
        logger.debug(EventType.MISC.getMarker(), "Setting fetch size {}", optionValue);
        Pattern pattern = Pattern.compile("(\\d+) MB");
        Matcher matcher = pattern.matcher(optionValue);
        if (matcher.find()) {
            chunkSize = Integer.parseInt(matcher.group(1)) * 1_000_000;
        } else {
            changedFetchSize = true;
            currentFetchSize = Integer.parseInt(optionValue);
        }
    }

    private void setInitFetchSize(String optionValue) {
        if (optionValue.isEmpty()) {
            return;
        }
        logger.debug(EventType.MISC.getMarker(), "Setting init fetch size {}", optionValue);
        initFetchSize = Integer.parseInt(optionValue);
    }
}

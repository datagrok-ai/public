package grok_connect.utils;

import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import grok_connect.GrokConnect;
import grok_connect.connectors_info.FuncCall;
import grok_connect.log.EventType;
import grok_connect.log.QueryLogger;
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
    private final Logger logger;
    private final QueryLogger<DataFrame> queryLogger;
    private final FuncCall query;
    private int currentFetchSize = MIN_FETCH_SIZE;
    private int initFetchSize = MIN_FETCH_SIZE;
    private int chunkSize = -1;
    private SchemeInfo schemeInfo;
    private ResultSet resultSet;
    private Connection connection;
    private boolean changedFetchSize;
    private boolean supportTransactions;
    private String initMessage;

    public QueryManager(String message, QueryLogger<DataFrame> queryLogger) {
        query = gson.fromJson(message, FuncCall.class);
        query.setParamValues();
        query.afterDeserialization();
        this.logger = queryLogger.getLogger();
        this.queryLogger = queryLogger;
        provider = GrokConnect.providerManager.getByName(query.func.connection.dataSource);
        initParams();
        if (dryRun)
            initMessage = message;
    }

    public void initResultSet(FuncCall query) throws ClassNotFoundException, GrokConnectException, QueryCancelledByUser, SQLException {
        logger.debug(EventType.CONNECTION_RECEIVE.getMarker(EventType.Stage.START), "Receiving connection to db");
        connection = provider.getConnection(query.func.connection);
        logger.debug(EventType.CONNECTION_RECEIVE.getMarker(EventType.Stage.END), "Connection was received");
        resultSet = provider.getResultSet(query, connection, logger, initFetchSize);
        supportTransactions = connection.getMetaData().supportsTransactions();
    }

    public void initScheme() throws QueryCancelledByUser, SQLException {
        if (resultSet == null) {
            logger.debug(EventType.SCHEME_INFO_INIT.getMarker(EventType.Stage.END), "resultSet is null, return empty schemeInfo");
            schemeInfo = new SchemeInfo(new ArrayList<>(), new ArrayList<>(), new ArrayList<>());
        } else {
            logger.debug(EventType.SCHEME_INFO_INIT.getMarker(EventType.Stage.START), "Initializing schemeInfo");
            schemeInfo = provider.resultSetScheme(query, resultSet, logger);
            logger.debug(EventType.SCHEME_INFO_INIT.getMarker(EventType.Stage.END), "Finished schemeInfo init");
        }
    }

    public void dryRun(boolean skipColumnFillingLog) throws QueryCancelledByUser, SQLException, IOException, GrokConnectException, ClassNotFoundException {
        // need to create new FuncCall on every run because state of params changes irrevocably after each run
        queryLogger.writeLog(false);
        FuncCall query = gson.fromJson(initMessage, FuncCall.class);
        query.setParamValues();
        initResultSet(query);
        initScheme();
        if (supportTransactions && changedFetchSize) {
            resultSet.setFetchSize(currentFetchSize);
        }
        queryLogger.writeLog(!skipColumnFillingLog);
        provider.getResultSetSubDf(query, resultSet, schemeInfo.columns,
                schemeInfo.supportedType, schemeInfo.initColumn, -1, logger, 1, true);
        closeConnection();
        queryLogger.writeLog(true);
    }

    public DataFrame getSubDF(int dfNumber) throws IOException, SQLException, QueryCancelledByUser {
        DataFrame df = new DataFrame();
        if (resultSet.isAfterLast() && !resultSet.isBeforeFirst())
            return df;

        schemeInfo.columns.forEach(Column::empty);

        if (!connection.isClosed() && !resultSet.isClosed()) {
            int rowsNumber = dfNumber == 1 ? initFetchSize : currentFetchSize;
            df =  provider.getResultSetSubDf(query, resultSet, schemeInfo.columns,
                    schemeInfo.supportedType, schemeInfo.initColumn, rowsNumber, logger, dfNumber, false);
        }

        if (dfNumber == 1 && supportTransactions && changedFetchSize) {
            logger.debug(EventType.MISC.getMarker(), "Manual fetch size was set for all chunks {}", currentFetchSize);
            resultSet.setFetchSize(currentFetchSize);
        } else if (!changedFetchSize)
            changeFetchSize(df, dfNumber);

        df.tags = new LinkedHashMap<>();
        df.tags.put(CHUNK_NUMBER_TAG, String.valueOf(dfNumber));
        return df;
    }

    public void closeConnection() throws SQLException {
        if (connection != null && !connection.isClosed()) {
            logger.debug(EventType.MISC.getMarker(), "Closing DB connection");
            if (!connection.getAutoCommit())
                connection.commit();
            provider.providerManager.getQueryMonitor().removeResultSet(query.id);
            connection.close();
            logger.debug(EventType.MISC.getMarker(), "DB connection was closed");
        } else
            provider.providerManager.getQueryMonitor().removeResultSet(query.id);
    }

    public boolean isResultSetInitialized() {
        return resultSet != null;
    }

    public FuncCall getQuery() {
        return query;
    }

    private void changeFetchSize(DataFrame df, int dfNumber) throws SQLException {
        if (supportTransactions && df.rowCount != 0) {
            currentFetchSize = getFetchSize(df);
            logger.debug(EventType.MISC.getMarker(dfNumber), "Calculated fetch size: {}", currentFetchSize);
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

    private void initParams() {
        Map<String, Object> params = new HashMap<>();
        // add meta params of query
        if (query.func.options != null)
            params.putAll(query.func.options);
        // aux params have preference over meta params
        if (query.func.aux != null) {
            for (String key: query.func.aux.keySet()) {
                Object value = query.func.aux.get(key);
                if (value != null && !value.toString().isEmpty())
                    params.put(key, value);
            }
        }
        setFetchSize(params.getOrDefault(FETCH_SIZE_KEY, "").toString());
        setInitFetchSize(params.getOrDefault(INIT_FETCH_SIZE_KEY, "").toString());
        dryRun = Boolean.parseBoolean(params.getOrDefault(DRY_RUN_KEY, false).toString());
    }

    private void setFetchSize(String optionValue) {
        if (optionValue == null || optionValue.isEmpty())
            return;
        Pattern pattern = Pattern.compile("(\\d+) MB");
        Matcher matcher = pattern.matcher(optionValue);
        if (matcher.find())
            chunkSize = Integer.parseInt(matcher.group(1)) * 1_000_000;
        else {
            changedFetchSize = true;
            currentFetchSize = Integer.parseInt(optionValue);
        }
    }

    private void setInitFetchSize(String optionValue) {
        if (optionValue == null || optionValue.isEmpty())
            return;
        logger.debug(EventType.MISC.getMarker(), "Setting init fetch size {}", optionValue);
        initFetchSize = Integer.parseInt(optionValue);
    }
}

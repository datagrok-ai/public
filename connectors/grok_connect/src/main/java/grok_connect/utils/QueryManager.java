package grok_connect.utils;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataProvider;
import grok_connect.connectors_info.FuncCall;
import grok_connect.handlers.QueryHandler;
import grok_connect.log.EventType;
import grok_connect.providers.JdbcDataProvider;
import grok_connect.resultset.ResultSetManager;
import grok_connect.table_mutation.MutationValidationException;
import grok_connect.table_mutation.TableMutation;
import grok_connect.table_query.TableQuery;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.slf4j.MDC;
import serialization.BigIntColumn;
import serialization.Column;
import serialization.DataFrame;
import serialization.Types;

public class QueryManager {
    public static final String CHUNK_NUMBER_TAG = "chunkNumber";
    public static final String FETCH_SIZE_KEY = "connectFetchSize";
    public static final String INIT_FETCH_SIZE_KEY = "initConnectFetchSize";
    public static final String INT8_AS_INT32_KEY = "int8AsInt32";
    public static final String COMPRESS_CHUNKS_KEY = "compressChunks";
    public static final String GZIP_LEVEL_KEY = "gzipLevel";
    public static final String COLUMN_TIMINGS_KEY = "columnTimings";
    public static final String ALLOW_COL_TYPE_CHANGE_TAG = "ALLOW_COL_TYPE_CHANGE";
    private static final int DEFAULT_GZIP_LEVEL = 1;
    private static final int MAX_CHUNK_SIZE_BYTES = 10_000_000;
    private static final int MAX_FETCH_SIZE = 500000;
    private static final int MIN_FETCH_SIZE = 100;
    private static final Logger LOGGER = LoggerFactory.getLogger(QueryManager.class);
    public boolean isDebug;
    public final JdbcDataProvider provider;
    private final FuncCall query;
    private final ResultSetManager resultSetManager;
    private int currentFetchSize = MIN_FETCH_SIZE;
    private int initFetchSize = MIN_FETCH_SIZE;
    private int chunkSize = -1;
    private ResultSet resultSet;
    private Connection connection;
    private boolean changedFetchSize;
    private String initMessage;
    private int columnCount;
    private boolean isFinished = false;
    private boolean supportsFetchSize = true;
    private boolean committed;
    private volatile float lastWireBytesPerRow;
    private boolean hasDowncast;

    public QueryManager(String message) {
        LOGGER.debug("Deserializing json call and preprocessing it...");
        query = GrokConnect.gson.fromJson(message, FuncCall.class);
        // mutations must never enter the query/dryRun path (dryRun executes twice and commits)
        if (query.func instanceof TableMutation)
            throw new MutationValidationException("Table mutations are not allowed on the query socket; use POST /mutate");
        query.setParamValues();
        query.afterDeserialization();
        isDebug = query.debugQuery;
        provider = GrokConnect.providerManager.getByName(query.func.connection.dataSource);
        resultSetManager = provider.getResultSetManager();
        initParams();
        if (isDebug)
            initMessage = message;
        LOGGER.debug("Deserialized and preprocessed call");
        processTableQuery(query);
    }

    private void processTableQuery(FuncCall query) {
        if (query.func instanceof TableQuery) {
            LOGGER.debug("Building table query...");
            query.func.query = provider.queryTableSql((TableQuery) query.func);
            String catalog = ((TableQuery) query.func).catalog;
            if (GrokConnectUtil.isNotEmpty(catalog))
                query.func.connection.parameters.put("db", catalog);
            LOGGER.debug("TableQuery was built");
        }
    }

    public void initResultSet(FuncCall query) throws GrokConnectException, QueryCancelledByUser, SQLException {
        committed = false; // dryRun re-inits the SAME instance: a stale flag would skip the real run's commit
        LOGGER.debug(EventType.CONNECTION_RECEIVE.getMarker(EventType.Stage.START), "Receiving connection to {} database...", provider.descriptor.type);
        connection = provider.getConnection(query.func.connection);
        LOGGER.debug(EventType.CONNECTION_RECEIVE.getMarker(EventType.Stage.END), "Received connection to {} database", provider.descriptor.type);
        resultSet = provider.getResultSet(query, connection, initFetchSize);
        if (resultSet == null) return;
        ResultSetMetaData metaData = resultSet.getMetaData();
        LOGGER.debug("Initializing ResultSet manager...");
        resultSetManager.init(metaData, currentFetchSize);
        LOGGER.debug("ResultSet manager was initialized");
        columnCount = metaData.getColumnCount();
    }

    public void dryRun(boolean skipColumnFillingLog) throws QueryCancelledByUser, SQLException, GrokConnectException {
        FuncCall query = GrokConnect.gson.fromJson(initMessage, FuncCall.class);
        query.setParamValues();
        processTableQuery(query);
        initResultSet(query);
        String sessionId = MDC.get(QueryHandler.CALL_ID_HEADER);
        if (skipColumnFillingLog)
            MDC.remove(QueryHandler.CALL_ID_HEADER);
        if (changedFetchSize)
            tryFetchSize(currentFetchSize);
        provider.getResultSetSubDf(query, resultSet, provider.getResultSetManager(), -1, columnCount,
                1, true);
        MDC.remove(QueryHandler.CALL_ID_HEADER);
        close(true);
        MDC.put(QueryHandler.CALL_ID_HEADER, sessionId);
    }

    public DataFrame getSubDF(int dfNumber) throws SQLException, QueryCancelledByUser {
        return getSubDF(dfNumber, 0);
    }

    /** {@code wireBytesPerRow} (the caller's snapshot of {@link #getWireBytesPerRow()}, 0 = unknown) sizes the chunk
     *  after this one, so the decision does not depend on which serialize has completed while this fetch runs. */
    public DataFrame getSubDF(int dfNumber, float wireBytesPerRow) throws SQLException, QueryCancelledByUser {
        DataFrame df = new DataFrame();
        if (!isFinished && !resultSet.isClosed() && !connection.isClosed()) {
            if (dfNumber != 1)
                resultSetManager.detach(currentFetchSize);
            int rowsNumber = dfNumber == 1 ? initFetchSize : currentFetchSize;
            df =  provider.getResultSetSubDf(query, resultSet, resultSetManager, rowsNumber, columnCount, dfNumber, false);
            if (df.rowCount == rowsNumber && supportsFetchSize) {
                if (dfNumber == 1 && changedFetchSize)
                    tryFetchSize(currentFetchSize);
                else if (!changedFetchSize)
                    changeFetchSize(df, wireBytesPerRow);
            }
            else {
                isFinished = true;
                LOGGER.info("Received all data");
            }
            tagChunk(df, dfNumber);
            if (dfNumber == 1)
                try {
                    SqlAnnotator.annotate(query.func, df);
                } catch (Exception ignore) {}
        }
        return df;
    }

    /** Commits once, right before the COMPLETED token is sent (WS frames are FIFO, so this
     *  happens-before Datlas resolving the caller's future); close(true) is a no-op afterwards. */
    public void commitPending() throws SQLException {
        if (connection != null && !connection.isClosed() && !connection.getAutoCommit() && !committed) {
            connection.commit();
            committed = true;
        }
    }

    public void close(boolean commit) throws SQLException {
        if (resultSet != null && !resultSet.isClosed())
            resultSet.close();
        if (connection != null && !connection.isClosed()) {
            LOGGER.debug("Closing DB connection...");
            if (!committed && !connection.getAutoCommit()) {
                if (commit)
                    connection.commit();
                else
                    provider.rollbackQuietly(connection); // failed rollback evicts — never pool a dirty connection
            }
            QueryMonitor.getInstance().removeResultSet(query.id);
            connection.close();
            LOGGER.debug("Closed DB connection");
        } else
            QueryMonitor.getInstance().removeResultSet(query.id);
    }

    public boolean isResultSetInitialized() {
        return resultSet != null;
    }

    public FuncCall getQuery() {
        return query;
    }

    /** True when the §6.2 post-hoc detector flagged this call (an update count was observed on a
     *  no-result-set statement while allowRawWrites auditing was requested — WO-B13). */
    public boolean isRawWriteDetected() {
        return Boolean.TRUE.equals(query.aux.get(DataProvider.RAW_WRITE_DETECTED));
    }

    public void reportSerialized(int bytes, int rows) {
        if (rows > 0)
            lastWireBytesPerRow = bytes / (float) rows;
    }

    public float getWireBytesPerRow() {
        return lastWireBytesPerRow;
    }

    void tagChunk(DataFrame df, int dfNumber) {
        df.setTag(CHUNK_NUMBER_TAG, String.valueOf(dfNumber));
        if (dfNumber == 1)
            hasDowncast = hasDowncastColumn(df);
        if (hasDowncast)
            df.setTag(ALLOW_COL_TYPE_CHANGE_TAG, "true");
    }

    static boolean int8AsInt32(Map<String, Object> options) {
        return String.valueOf(options.get(INT8_AS_INT32_KEY)).equals("true");
    }

    public static boolean compressChunks(Map<String, Object> options) {
        return String.valueOf(options.get(COMPRESS_CHUNKS_KEY)).equals("true");
    }

    public static boolean columnTimings(Map<String, Object> options) {
        return String.valueOf(options.get(COLUMN_TIMINGS_KEY)).equals("true");
    }

    // Mirrors Datlas (external_provider.dart): chunk 1 is gzipped only when initConnectFetchSize
    // was given and is not the number 100 (the string "100" counts as given).
    public static boolean gzipChunk(Map<String, Object> options, int dfNumber) {
        if (!compressChunks(options))
            return false;
        Object initFetchSize = options.get(INIT_FETCH_SIZE_KEY);
        return dfNumber > 1 || (initFetchSize != null
                && !(initFetchSize instanceof Number && ((Number) initFetchSize).doubleValue() == 100));
    }

    public static int gzipLevel(Map<String, Object> options) {
        Object level = options.get(GZIP_LEVEL_KEY);
        if (level == null)
            return DEFAULT_GZIP_LEVEL;
        int value = level instanceof Number ? ((Number) level).intValue() : Double.valueOf(level.toString()).intValue();
        return Math.min(Math.max(1, value), 9);
    }

    private static boolean hasDowncastColumn(DataFrame df) {
        for (Column<?> column : df.getColumns())
            if (column instanceof BigIntColumn && column.getType().equals(Types.INT))
                return true;
        return false;
    }

    private void changeFetchSize(DataFrame df, float wireBytesPerRow) {
        if (!provider.descriptor.type.equals("Virtuoso")) {
            LOGGER.debug("Calculating dynamically next fetch size...");
            currentFetchSize = getFetchSize(df, wireBytesPerRow);
            tryFetchSize(currentFetchSize);
        }
    }

    int getFetchSize(DataFrame df, float wireBytesPerRow) {
        int fetchSize = 0;
        if (df.rowCount != 0) {
            int maxChunkSize = chunkSize != -1 ? chunkSize : MAX_CHUNK_SIZE_BYTES;
            float bytesPerRow = wireBytesPerRow > 0 ? wireBytesPerRow : (float) (df.memoryInBytes() / df.rowCount);
            fetchSize = Math.round(maxChunkSize / bytesPerRow);
        }
        return Math.min(Math.max(MIN_FETCH_SIZE, fetchSize), MAX_FETCH_SIZE);
    }

    private void initParams() {
        query.options.entrySet().removeIf(e -> e.getValue() == null);
        setFetchSize(query.options.getOrDefault(FETCH_SIZE_KEY, "").toString());
        setInitFetchSize(query.options.getOrDefault(INIT_FETCH_SIZE_KEY, "").toString());
        resultSetManager.setInt8AsInt32(int8AsInt32(query.options));
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
            currentFetchSize = Double.valueOf(optionValue).intValue();
            LOGGER.info("Fetch size was set to {} for all chunks except initial", currentFetchSize);
        }
    }

    private void setInitFetchSize(String optionValue) {
        if (optionValue == null || optionValue.isEmpty()) {
            LOGGER.info("Default init fetch size of {} will be used", initFetchSize);
            return;
        }
        try {
            initFetchSize = Double.valueOf(optionValue).intValue();
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException(INIT_FETCH_SIZE_KEY + " must be a row count, got '" + optionValue
                    + "' (the \"N MB\" form applies to " + FETCH_SIZE_KEY + " only)");
        }
        LOGGER.info("Init fetch size was set to {}", initFetchSize);
    }

    private void tryFetchSize(int fetchSize) {
        try {
            resultSet.setFetchSize(fetchSize);
            LOGGER.info("Fetch size was set to {}", currentFetchSize);
        } catch (SQLException e) {
            LOGGER.info("Provider {} doesn't support fetch size change", provider.descriptor.type);
            supportsFetchSize = false;
        }
    }
}

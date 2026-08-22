package grok_connect.handlers;

import grok_connect.connectors_info.DataQueryRunResult;
import grok_connect.log.EventType;
import grok_connect.providers.DatabricksProvider;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.QueryChunkNotSent;
import grok_connect.utils.QueryManager;
import grok_connect.utils.QueryMonitor;
import grok_connect.table_mutation.MutationManager;
import grok_connect.table_mutation.MutationResult;
import grok_connect.table_mutation.MutationValidationException;
import org.eclipse.jetty.websocket.api.Session;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.*;
import java.util.stream.Collectors;
import java.util.zip.GZIPOutputStream;
import grok_connect.GrokConnect;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.slf4j.MDC;
import org.slf4j.Marker;
import serialization.BufferAccessor;
import serialization.Column;
import serialization.DataFrame;
import serialization.TablesBlob;

public class SessionHandler {
    public static final Logger LOGGER = LoggerFactory.getLogger(SessionHandler.class);
    private static final List<String> DEFAULT_LOG_LEVELS = Arrays.asList("info", "warn", "error");
    private static final String COMPLETED_OK = "COMPLETED_OK";
    private static final String MESSAGE_START = "QUERY";
    private static final String OK_RESPONSE = "PART OK";
    private static final String END_MESSAGE = "EOF";
    private static final String SIZE_RECEIVED_MESSAGE = "DATAFRAME PART SIZE RECEIVED";
    private static final String MUTATION_START = "MUTATION ";
    private static final String MUTATION_EOF = "MUTATION EOF";
    private static final String MUTATION_READY = "MUTATION READY";
    private static final long MUTATION_IDLE_TIMEOUT_MS = 10 * 60 * 1000;
    private static final String RESULT_PREFIX = "RESULT ";
    private static final long CANCEL_JOIN_TIMEOUT_MS = 5000;
    static boolean PIPELINE = Boolean.parseBoolean(System.getProperty("grok.connect.pipelineFetch", "true"));
    private final Session session;
    private final boolean skipLogging;
    private CompletableFuture<Frame> completableFuture;
    private CompletableFuture<DataFrame> fetchFuture;
    private Frame frame;
    private Boolean firstTry = true;
    private Boolean oneDfSent = false;
    private int dfNumber = 1;
    private QueryManager queryManager;
    private MutationManager mutationManager;
    private boolean completedOk = false;

    SessionHandler(Session session, boolean skipLog) {
        this.session = session;
        this.skipLogging = skipLog;
    }

    public void onError(Throwable err) {
        if (err instanceof OutOfMemoryError) {
            // guess it won't work because there is no memory left!
            GrokConnect.needToReboot = true;
        }
        if (err.getClass().equals(GrokConnectException.class) && err.getCause() != null)
            err = err.getCause();
        String message = err.getMessage();
        // todo: improve this
        if (queryManager != null && queryManager.provider.getClass().equals(DatabricksProvider.class))
            message = DatabricksProvider.simplifyDatabricksError(message);
        String stackTrace = Arrays.stream(err.getStackTrace()).map(StackTraceElement::toString)
                .collect(Collectors.joining(System.lineSeparator()));
        LOGGER.error(EventType.ERROR.getMarker(), message, err);
        DataQueryRunResult result = new DataQueryRunResult();
        result.errorMessage = message;
        result.errorStackTrace = stackTrace;
        session.getRemote().sendStringByFuture(String.format("ERROR: %s", GrokConnect.gson.toJson(result)));
        closeQueryManagerQuietly();
        session.close();
    }

    public void onMessage(String message) throws Throwable {
        try {
            onMessageInternal(message);
        } catch (Throwable t) {
            // Ensure the DB connection is released even if onClose never fires
            // (half-open WS, datlas crash mid-stream, errors before COMPLETED_OK).
            closeQueryManagerQuietly();
            throw t;
        }
    }

    public void onBinary(byte[] payload, int offset, int len) throws Throwable {
        try {
            if (mutationManager == null)
                throw new MutationValidationException("Received a binary frame outside a mutation session");
            // Jetty reuses the frame buffer after this call returns — copy before handing to the loader.
            byte[] chunk = Arrays.copyOfRange(payload, offset, offset + len);
            mutationManager.feed(chunk);
            session.getRemote().sendStringByFuture(OK_RESPONSE);
        } catch (Throwable t) {
            closeQueryManagerQuietly();
            throw t;
        }
    }

    private void onMessageInternal(String message) throws Throwable {
        if (message.equals(MUTATION_EOF)) {
            if (mutationManager == null)
                throw new MutationValidationException("No mutation session to finish");
            MutationResult result = mutationManager.finish();
            session.getRemote().sendStringByFuture(RESULT_PREFIX + GrokConnect.gson.toJson(result));
            session.close();
            return;
        }
        if (message.startsWith(MUTATION_START)) {
            if (queryManager != null)
                throw new MutationValidationException("Query session cannot process mutation messages");
            if (mutationManager != null)
                throw new MutationValidationException("Mutation already started for this session");
            LOGGER.debug("Received bulk mutation header from the server");
            // Queries keep the connect-time idle timeout of 0 (a long-running statement legitimately
            // produces no WS traffic), but a mutation session receives chunks continuously — prolonged
            // silence means the sender died. Without a timeout a half-open socket would hold the open
            // transaction (and its pooled connection) forever; on timeout Jetty fires onError/onClose,
            // which reach mutationManager.abort() -> rollback before the connection returns to the pool.
            session.setIdleTimeout(MUTATION_IDLE_TIMEOUT_MS);
            mutationManager = new MutationManager(message.substring(MUTATION_START.length()));
            mutationManager.start();
            session.getRemote().sendStringByFuture(MUTATION_READY);
            return;
        }
        if (message.startsWith(MESSAGE_START)) {
            if (mutationManager != null)
                throw new MutationValidationException("Mutation session cannot process query messages");
            LOGGER.debug("Received message with json call from the server");
            message = message.substring(6);
            queryManager = new QueryManager(message);
            if (queryManager.isDebug) {
                LOGGER.info(EventType.DRY_RUN.getMarker(EventType.Stage.START), "Running dry run...");
                for (int i = 0; i < 2; i++) {
                    LOGGER.debug("Running #{} dry run {} logging dataframe filling times...", i + 1, i == 0 ? "without" : "with");
                    queryManager.dryRun(i == 0);
                    LOGGER.debug("Finished #{} dry run.", i + 1);
                }
                LOGGER.info(EventType.DRY_RUN.getMarker(EventType.Stage.END), "Dry run finished");
            }

            queryManager.initResultSet(queryManager.getQuery());
            frame = serialize(queryManager.isResultSetInitialized() ? queryManager.getSubDF(dfNumber) : new DataFrame(), dfNumber);
            if (PIPELINE && queryManager.isResultSetInitialized()) {
                float bytesPerRow = queryManager.getWireBytesPerRow();
                fetchFuture = async(MDC.getCopyOfContextMap(), () -> queryManager.getSubDF(2, bytesPerRow));
            }
        }
        else if (message.startsWith(SIZE_RECEIVED_MESSAGE)) {
            // A re-announce (firstTry retry) must not submit a second read of the same chunk.
            if (completableFuture == null && queryManager.isResultSetInitialized()) {
                int next = dfNumber + 1;
                Map<String, String> context = MDC.getCopyOfContextMap();
                // Snapshot of chunk dfNumber's serialized bytes/row: the fetch scheduled here sizes its successor from
                // it, not from whatever serialize(next) has written by the time that fetch ends.
                float bytesPerRow = queryManager.getWireBytesPerRow();
                if (fetchFuture == null)
                    completableFuture = async(context, () -> toFrame(queryManager.getSubDF(next, bytesPerRow), next));
                else {
                    // serialize(next) and fetch(next + 1) run on separate data (detach); fetch(next + 2) is chained only
                    // by the next SIZE RECEIVED, which bounds live data to two DataFrames and two Frames.
                    CompletableFuture<DataFrame> fetched = fetchFuture;
                    completableFuture = fetched.thenCompose(df -> async(context, () -> toFrame(df, next)));
                    fetchFuture = fetched.thenCompose(df -> df.rowCount == 0 ? fetched : async(context, () -> queryManager.getSubDF(next + 1, bytesPerRow)));
                }
            }
            Marker start = EventType.SOCKET_BINARY_DATA_EXCHANGE.getMarker(dfNumber, EventType.Stage.START);
            LOGGER.debug(start, "Sending binary data with id {} to the server...", dfNumber);
            session.getRemote().sendBytesByFuture(ByteBuffer.wrap(frame.bytes));
            return;
        }
        else if (message.startsWith(COMPLETED_OK)) {
            completedOk = true;
            session.getRemote().sendStringByFuture(END_MESSAGE);
            session.close();
            return;
        }
        else {
            if (!message.equals(OK_RESPONSE)) {
                if (!firstTry)
                    throw new QueryChunkNotSent();
                else
                    firstTry = false;
            }
            else {
                firstTry = true;
                oneDfSent = true;
                if (queryManager.isResultSetInitialized()) {
                    LOGGER.debug("-- GrokConnect thread is doing nothing. Blocking to receive next dataframe --");
                    try {
                        frame = completableFuture.join();
                    } catch (CompletionException e) {
                        // rethrow exact cause
                        throw e.getCause();
                    }
                    completableFuture = null;
                    LOGGER.debug("-- Received next dataframe from executing thread --");
                    if (frame.rowCount != 0) dfNumber++;
                }
            }
        }
        if (frame != null && (frame.rowCount != 0 || !oneDfSent)) {
            LOGGER.debug(EventType.CHECKSUM_SEND.getMarker(dfNumber, EventType.Stage.START), "Sending checksum message. Data size: {}", frame.bytes.length);
            session.getRemote().sendStringByFuture(String.format("DATAFRAME PART SIZE: %d%s", frame.bytes.length, frame.gzipped ? " gzip=true" : ""));
        }
        else {
            // Commit BEFORE the COMPLETED token: Datlas resolves the caller's future on it, and WS
            // frames are FIFO — commit-before-send is the happens-before edge that closes the
            // commit-visibility race. Safe only here: the state machine reaches this branch with a
            // fully-drained resultset. A failed commit propagates as a query ERROR instead of being
            // swallowed post-success in onClose.
            if (queryManager != null)
                queryManager.commitPending();
            // The RAW_WRITE token carries the §6.2 post-hoc raw-write detection (WO-B13) to Datlas,
            // which audit-logs it; older servers only read the chunk-count token and ignore the rest.
            String rawWrite = queryManager != null && queryManager.isRawWriteDetected() ? " RAW_WRITE" : "";
            session.getRemote().sendStringByFuture(String.format("COMPLETED %s%s", dfNumber, rawWrite));
        }
    }

    public void onClose() throws SQLException {
        cancelTasks();
        if (queryManager != null)
            queryManager.close(completedOk);
        if (mutationManager != null)
            mutationManager.abort(); // no-op if finish() already committed and closed
    }

    private void cancelTasks() {
        CompletableFuture<?> fetching = fetchFuture != null ? fetchFuture : completableFuture;
        if (fetching != null && queryManager != null) {
            // cancel(true) never interrupts a running pool task: flag the query so the row loop throws
            // QueryCancelledByUser, and wait for that before the result set and connection are closed under it.
            QueryMonitor.getInstance().addCancelledResultSet(queryManager.getQuery().id);
            try {
                fetching.get(CANCEL_JOIN_TIMEOUT_MS, TimeUnit.MILLISECONDS);
            } catch (Exception ignore) {}
        }
        if (completableFuture != null)
            completableFuture.cancel(true);
        if (fetchFuture != null)
            fetchFuture.cancel(true);
    }

    private static <T> CompletableFuture<T> async(Map<String, String> context, Callable<T> body) {
        return GrokConnect.submitPoolTask(() -> {
            try {
                if (context != null)
                    MDC.setContextMap(context);
                return body.call();
            } catch (Exception e) {
                throw new CompletionException(e);
            } finally {
                MDC.clear();
            }
        });
    }

    private Frame toFrame(DataFrame df, int chunkNumber) throws IOException {
        return df.rowCount != 0 ? serialize(df, chunkNumber) : new Frame(null, false, 0);
    }

    private void closeQueryManagerQuietly() {
        cancelTasks();
        if (queryManager != null) {
            try {
                queryManager.close(false);
            } catch (Throwable t) {
                LOGGER.warn("Failed to close QueryManager during cleanup", t);
            }
        }
        if (mutationManager != null) {
            try {
                mutationManager.abort();
            } catch (Throwable t) {
                LOGGER.warn("Failed to abort MutationManager during cleanup", t);
            }
        }
    }

    Frame serialize(DataFrame df, int chunkNumber) throws IOException {
        LOGGER.debug(EventType.DATAFRAME_TO_BYTEARRAY_CONVERSION.getMarker(chunkNumber, EventType.Stage.START), "Converting DataFrame with id {} to binary data...", chunkNumber);
        TablesBlob blob = df.toBlob();
        byte[] bytes = blob.toByteArray();
        queryManager.reportSerialized(bytes.length, df.rowCount);
        LOGGER.debug(EventType.DATAFRAME_TO_BYTEARRAY_CONVERSION.getMarker(chunkNumber, EventType.Stage.END), "Converted DataFrame with id {} to binary data", chunkNumber);
        Map<String, Object> options = queryManager.getQuery().options;
        if (queryManager.isDebug)
            LOGGER.debug(EventType.MISC.getMarker(chunkNumber), "COLUMN SIZES {}", GrokConnect.gson.toJson(columnSizes(df, blob, bytes, QueryManager.columnTimings(options))));
        if (!QueryManager.gzipChunk(options, chunkNumber))
            return new Frame(bytes, false, df.rowCount);
        LOGGER.debug(EventType.COMPRESSION.getMarker(chunkNumber, EventType.Stage.START), "Compressing binary data with id {}...", chunkNumber);
        byte[] gzipped = gzip(bytes, QueryManager.gzipLevel(options));
        LOGGER.debug(EventType.COMPRESSION.getMarker(chunkNumber, EventType.Stage.END), "Compressed binary data with id {}: {} -> {} bytes", chunkNumber, bytes.length, gzipped.length);
        return new Frame(gzipped, true, df.rowCount);
    }

    static byte[] gzip(byte[] data, int level) throws IOException {
        ByteArrayOutputStream out = new ByteArrayOutputStream(data.length / 3);
        try (GZIPOutputStream gz = new GZIPOutputStream(out, 65536) {{ def.setLevel(level); }}) {
            gz.write(data);
        }
        return out.toByteArray();
    }

    static final class Frame {
        final byte[] bytes;
        final boolean gzipped;
        final int rowCount;

        Frame(byte[] bytes, boolean gzipped, int rowCount) {
            this.bytes = bytes;
            this.gzipped = gzipped;
            this.rowCount = rowCount;
        }
    }

    private static List<Map<String, Object>> columnSizes(DataFrame df, TablesBlob blob, byte[] bytes, boolean timings) throws IOException {
        int[] offsets = blob.getColumnsOffsets(0);
        BufferAccessor reader = new BufferAccessor(bytes);
        List<Map<String, Object>> result = new ArrayList<>();
        for (int i = 0; i < offsets.length; i++) {
            int end = i + 1 < offsets.length ? offsets[i + 1] : blob.getLength();
            reader.bufPos = offsets[i];
            ByteArrayOutputStream gz = new ByteArrayOutputStream();
            try (GZIPOutputStream out = new GZIPOutputStream(gz)) {
                out.write(bytes, offsets[i], end - offsets[i]);
            }
            Column<?> col = df.getColumn(i);
            Map<String, Object> column = new LinkedHashMap<>();
            column.put("name", col.getName());
            column.put("type", col.getType());
            column.put("enc", reader.readColumnEncoderId());
            column.put("bytes", end - offsets[i]);
            column.put("gz", gz.size());
            if (timings) {
                // Opt-in re-encode of the column alone: the serialization writer exposes no per-column timing.
                long start = System.nanoTime();
                col.encode(new BufferAccessor());
                column.put("ms", Math.round((System.nanoTime() - start) / 1e5) / 10.0);
            }
            result.add(column);
        }
        return result;
    }

    public boolean skipLog(String level) {
        return skipLogging && !DEFAULT_LOG_LEVELS.contains(level);
    }

    public Session getSession() {
        return this.session;
    }
}

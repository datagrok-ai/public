package grok_connect.handlers;

import grok_connect.connectors_info.DataQueryRunResult;
import grok_connect.log.EventType;
import grok_connect.providers.DatabricksProvider;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.QueryCancelledByUser;
import grok_connect.utils.QueryChunkNotSent;
import grok_connect.utils.QueryManager;
import grok_connect.table_mutation.MutationManager;
import grok_connect.table_mutation.MutationResult;
import grok_connect.table_mutation.MutationValidationException;
import org.eclipse.jetty.websocket.api.Session;
import java.nio.ByteBuffer;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.*;
import java.util.stream.Collectors;
import grok_connect.GrokConnect;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.slf4j.MDC;
import org.slf4j.Marker;
import serialization.DataFrame;

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
    private static final String RESULT_PREFIX = "RESULT ";
    private final Session session;
    private final boolean skipLogging;
    private CompletableFuture<DataFrame> completableFuture;
    private DataFrame dataFrame;
    private Boolean firstTry = true;
    private Boolean oneDfSent = false;
    private int dfNumber = 1;
    private byte[] bytes;
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
            dataFrame = queryManager.isResultSetInitialized() ? queryManager.getSubDF(dfNumber) : new DataFrame();
        }
        else if (message.startsWith(SIZE_RECEIVED_MESSAGE)) {
            Map<String, String> copyOfContextMap = MDC.getCopyOfContextMap();
            completableFuture = GrokConnect.submitPoolTask(() -> {
                try {
                    if (copyOfContextMap != null)
                        MDC.setContextMap(copyOfContextMap);
                    return queryManager.getSubDF(dfNumber + 1);
                } catch (SQLException | QueryCancelledByUser e) {
                    throw new CompletionException(e);
                } finally {
                    MDC.clear();
                }
            });
            Marker start = EventType.SOCKET_BINARY_DATA_EXCHANGE.getMarker(dfNumber, EventType.Stage.START);
            LOGGER.debug(start, "Sending binary data with id {} to the server...", dfNumber);
            session.getRemote().sendBytesByFuture(ByteBuffer.wrap(bytes));
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
                        dataFrame = completableFuture.join();
                    } catch (CompletionException e) {
                        // rethrow exact cause
                        throw e.getCause();
                    }
                    LOGGER.debug("-- Received next dataframe from executing thread --");
                    if (dataFrame.rowCount != 0) dfNumber++;
                }
            }
        }
        if (dataFrame != null && (dataFrame.rowCount != 0 || !oneDfSent)) {
            Marker start = EventType.DATAFRAME_TO_BYTEARRAY_CONVERSION.getMarker(dfNumber, EventType.Stage.START);
            Marker finish = EventType.DATAFRAME_TO_BYTEARRAY_CONVERSION.getMarker(dfNumber, EventType.Stage.END);
            LOGGER.debug(start, "Converting DataFrame with id {} to binary data...", dfNumber);
            bytes = dataFrame.toByteArray();
            dataFrame = null;
            LOGGER.debug(finish, "Converted DataFrame with id {} to binary data", dfNumber);
            LOGGER.debug(EventType.CHECKSUM_SEND.getMarker(dfNumber, EventType.Stage.START), "Sending checksum message. Data size: {}", bytes.length);
            session.getRemote().sendStringByFuture(String.format("DATAFRAME PART SIZE: %d", bytes.length));
        }
        else {
            session.getRemote().sendStringByFuture(String.format("COMPLETED %s", dfNumber));
        }
    }

    public void onClose() throws SQLException {
        if (queryManager != null)
            queryManager.close(completedOk);
        if (mutationManager != null)
            mutationManager.abort(); // no-op if finish() already committed and closed
    }

    private void closeQueryManagerQuietly() {
        if (completableFuture != null && !completableFuture.isDone()) {
            completableFuture.cancel(true);
        }
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

    public boolean skipLog(String level) {
        return skipLogging && !DEFAULT_LOG_LEVELS.contains(level);
    }

    public Session getSession() {
        return this.session;
    }
}

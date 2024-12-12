package grok_connect.handlers;

import grok_connect.connectors_info.DataQueryRunResult;
import grok_connect.log.EventType;
import grok_connect.log.QueryLogger;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.QueryCancelledByUser;
import grok_connect.utils.QueryChunkNotSent;
import grok_connect.utils.QueryManager;
import org.eclipse.jetty.websocket.api.Session;
import java.nio.ByteBuffer;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.concurrent.*;
import java.util.stream.Collectors;
import grok_connect.GrokConnect;
import org.slf4j.Logger;
import org.slf4j.Marker;
import serialization.DataFrame;

public class SessionHandler {
    private static final String COMPLETED_OK = "COMPLETED_OK";
    private static final String MESSAGE_START = "QUERY";
    private static final String OK_RESPONSE = "DATAFRAME PART OK";
    private static final String END_MESSAGE = "EOF";
    private static final String SIZE_RECEIVED_MESSAGE = "DATAFRAME PART SIZE RECEIVED";
    private final Session session;
    private CompletableFuture<DataFrame> completableFuture;
    private final QueryLogger queryLogger;
    private final Logger logger;
    private DataFrame dataFrame;
    private Boolean firstTry = true;
    private Boolean oneDfSent = false;
    private int dfNumber = 1;
    private byte[] bytes;
    private QueryManager queryManager;

    SessionHandler(Session session, QueryLogger queryLogger) {
        this.queryLogger = queryLogger;
        this.logger = queryLogger.getLogger();
        this.session = session;
        logger.info("GrokConnect version: {}", GrokConnect.properties.getProperty("version"));
    }

    public void onError(Throwable err) {
        if (err instanceof OutOfMemoryError) {
            // guess it won't work because there is no memory left!
            GrokConnect.needToReboot = true;
        }
        if (err.getClass().equals(GrokConnectException.class))
            err = err.getCause();
        String message = err.getMessage();
        String stackTrace = Arrays.stream(err.getStackTrace()).map(StackTraceElement::toString)
                .collect(Collectors.joining(System.lineSeparator()));
        logger.error(EventType.ERROR.getMarker(), message, err);
        DataQueryRunResult result = new DataQueryRunResult();
        result.errorMessage = message;
        result.errorStackTrace = stackTrace;
        session.getRemote().sendStringByFuture(String.format("ERROR: %s", GrokConnect.gson.toJson(result)));
        session.close();
    }

    public void onMessage(String message) throws Throwable {
        if (message.startsWith(MESSAGE_START)) {
            logger.debug("Received message with json call from the server");
            message = message.substring(6);
            queryManager = new QueryManager(message, queryLogger);
            if (queryManager.isDebug) {
                logger.info(EventType.DRY_RUN.getMarker(EventType.Stage.START), "Running dry run...");
                for (int i = 0; i < 2; i++) {
                    logger.debug("Running #{} dry run {} logging dataframe filling times...", i + 1, i == 0 ? "without" : "with");
                    queryManager.dryRun(i == 0);
                    logger.debug("Finished #{} dry run.", i + 1);
                }
                logger.info(EventType.DRY_RUN.getMarker(EventType.Stage.END), "Dry run finished");
            }

            queryManager.initResultSet(queryManager.getQuery());
            dataFrame = queryManager.isResultSetInitialized() ? queryManager.getSubDF(dfNumber) : new DataFrame();
        }
        else if (message.startsWith(SIZE_RECEIVED_MESSAGE)) {
            completableFuture = GrokConnect.submitPoolTask(() -> {
                try {
                    return queryManager.getSubDF(dfNumber + 1);
                } catch (SQLException | QueryCancelledByUser e) {
                    throw new CompletionException(e);
                }
            });
            Marker start = EventType.SOCKET_BINARY_DATA_EXCHANGE.getMarker(dfNumber, EventType.Stage.START);
            logger.debug(start, "Sending binary data with id {} to the server...", dfNumber);
            session.getRemote().sendBytesByFuture(ByteBuffer.wrap(bytes));
            return;
        }
        else if (message.startsWith(COMPLETED_OK)) {
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
                    logger.debug("-- GrokConnect thread is doing nothing. Blocking to receive next dataframe --");
                    try {
                        dataFrame = completableFuture.join();
                    } catch (CompletionException e) {
                        // rethrow exact cause
                        throw e.getCause();
                    }
                    logger.debug("-- Received next dataframe from executing thread --");
                    if (dataFrame.rowCount != 0) dfNumber++;
                }
            }
        }
        if (dataFrame != null && (dataFrame.rowCount != 0 || !oneDfSent)) {
            Marker start = EventType.DATAFRAME_TO_BYTEARRAY_CONVERSION.getMarker(dfNumber, EventType.Stage.START);
            Marker finish = EventType.DATAFRAME_TO_BYTEARRAY_CONVERSION.getMarker(dfNumber, EventType.Stage.END);
            logger.debug(start, "Converting DataFrame with id {} to binary data...", dfNumber);
            bytes = dataFrame.toByteArray();
            dataFrame = null;
            logger.debug(finish, "Converted DataFrame with id {} to binary data", dfNumber);
            logger.debug(EventType.CHECKSUM_SEND.getMarker(dfNumber, EventType.Stage.START), "Sending checksum message. Data size: {}", bytes.length);
            session.getRemote().sendStringByFuture(String.format("DATAFRAME PART SIZE: %d", bytes.length));
        }
        else {
            session.getRemote().sendStringByFuture(String.format("COMPLETED %s", dfNumber));
        }
    }

    public void onClose() throws SQLException {
        if (queryManager != null) queryManager.close();
        queryLogger.closeLogger();
    }
}

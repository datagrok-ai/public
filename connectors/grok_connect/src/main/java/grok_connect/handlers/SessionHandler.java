package grok_connect.handlers;

import com.google.gson.Gson;
import grok_connect.connectors_info.DataQueryRunResult;
import grok_connect.connectors_info.FuncCall;
import grok_connect.log.EventType;
import grok_connect.log.QueryLogger;
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
    private static final String SIZE_RECIEVED_MESSAGE = "DATAFRAME PART SIZE RECEIVED";
    private final Session session;
    private final ExecutorService threadPool;
    private final QueryLogger<DataFrame> queryLogger;
    private final Logger logger;
    private Future<DataFrame> fdf;
    private DataFrame dataFrame;
    private Boolean firstTry = true;
    private Boolean oneDfSent = false;
    private int dfNumber = 1;
    private byte[] bytes;
    private QueryManager queryManager;

    SessionHandler(Session session, QueryLogger<DataFrame> queryLogger) {
        threadPool = Executors.newCachedThreadPool();
        this.queryLogger = queryLogger;
        this.logger = queryLogger.getLogger();
        this.session = session;
        logger.info(EventType.MISC.getMarker(), "GrokConnect version: {}", GrokConnect.properties.getProperty("version"));
    }

    public void onError(Throwable err) {
        Throwable cause = err.getCause() == null ? err : err.getCause(); // we need to get cause, because error is wrapped by runtime exception
        if (cause instanceof OutOfMemoryError) {
            // guess it won't work because there is no memory left!
            GrokConnect.needToReboot = true;
        }
        String message = cause.getMessage();
        String stackTrace = Arrays.stream(cause.getStackTrace()).map(StackTraceElement::toString)
                .collect(Collectors.joining(System.lineSeparator()));
        logger.error(EventType.ERROR.getMarker(), cause.toString(), cause);
        logger.error(EventType.ERROR.getMarker(), String.format("STACK_TRACE: %s", stackTrace));
        DataQueryRunResult result = new DataQueryRunResult();
        result.errorMessage = message;
        result.errorStackTrace = stackTrace;
        session.getRemote().sendStringByFuture(String.format("ERROR: %s", new Gson().toJson(result)));
        sendLog();
        session.close();
    }

    public void onMessage(String message) throws Throwable {
        logger.debug(EventType.MISC.getMarker(), "Received message {}", message.startsWith("QUERY") ? "QUERY" : message);
        if (message.startsWith(MESSAGE_START)) {
            message = message.substring(6);
            queryManager = new QueryManager(message, queryLogger);
            if (queryManager.dryRun) {
                logger.info(EventType.DRY_RUN.getMarker(EventType.Stage.START), "Running dry run");
                for (int i = 0; i < 2; i++)
                    queryManager.dryRun(i == 0);
                logger.info(EventType.DRY_RUN.getMarker(EventType.Stage.END), "Dry run finished");
            }
            queryManager.initResultSet(queryManager.getQuery());
            if (queryManager.isResultSetInitialized()) {
                dataFrame = queryManager.getSubDF(dfNumber);
            } else {
                dataFrame = new DataFrame();
            }
        } else if (message.startsWith(SIZE_RECIEVED_MESSAGE)) {
            fdf = threadPool.submit(() -> queryManager.getSubDF(dfNumber + 1));
            Marker start = EventType.DATA_SEND.getMarker(dfNumber, EventType.Stage.START);
            logger.debug(start, "Sending binary dataframe with id {}", dfNumber);
            session.getRemote().sendBytesByFuture(ByteBuffer.wrap(bytes));
            return;
        } else if (message.startsWith(COMPLETED_OK)) {
            sendLog();
            session.getRemote().sendStringByFuture(END_MESSAGE);
            session.close();
            return;
        } else {
            if (!message.equals(OK_RESPONSE)) {
                if (!firstTry)
                    throw new QueryChunkNotSent();
                else {
                    firstTry = false;
                }
            }
            else {
                firstTry = true;
                oneDfSent = true;
                if (queryManager.isResultSetInitialized()) {
                    dataFrame = fdf.get();
                    if (dataFrame != null && dataFrame.rowCount != 0) {
                        dfNumber++;
                    }
                }
            }
        }
        if (dataFrame != null && (dataFrame.rowCount != 0 || !oneDfSent)) {
            Marker start = EventType.DATAFRAME_TO_BYTEARRAY_CONVERSION.getMarker(dfNumber, EventType.Stage.START);
            Marker finish = EventType.DATAFRAME_TO_BYTEARRAY_CONVERSION.getMarker(dfNumber, EventType.Stage.END);
            logger.debug(start, "Converting dataframe to byteArray");
            bytes = dataFrame.toByteArray();
            logger.debug(finish, "DataFrame with id {} was converted to byteArray", dfNumber);
            logger.debug(EventType.CHECKSUM_SEND.getMarker(dfNumber, EventType.Stage.START), "Sending checkSum message with bytes length of {}", bytes.length);
            session.getRemote().sendStringByFuture(String.format("DATAFRAME PART SIZE: %d", bytes.length));
        } else {
            session.getRemote().sendStringByFuture(String.format("COMPLETED %s", dfNumber));
        }
    }

    public void onClose() throws SQLException {
        threadPool.shutdown();
        queryLogger.closeLogger();
        queryManager.closeConnection();
    }

    private void sendLog() {
        FuncCall query = queryManager.getQuery();
         if (query.debugQuery) {
            logger.debug(EventType.LOG_SEND.getMarker(EventType.Stage.START), "Converting logs to binary data and sending them to the server");
            DataFrame logs = queryLogger.dumpLogMessages();
            session.getRemote().sendBytesByFuture(ByteBuffer.wrap(logs.toByteArray()));
        }
    }
}

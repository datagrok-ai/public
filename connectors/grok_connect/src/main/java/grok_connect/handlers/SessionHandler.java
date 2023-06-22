package grok_connect.handlers;

import ch.qos.logback.classic.Level;
import com.google.gson.Gson;
import grok_connect.connectors_info.DataQueryRunResult;
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
    private int dryRunDiff = 0;
    private int dryRunCycles = 2;
    private byte[] bytes;
    private QueryManager queryManager;

    SessionHandler(Session session, QueryLogger<DataFrame> queryLogger) {
        threadPool = Executors.newCachedThreadPool();
        this.queryLogger = queryLogger;
        this.logger = queryLogger.getLogger();
        this.session = session;
        logger.info(EventType.MISC.getMarker(), "GrokConnect version: {}", GrokConnect.properties.getProperty("version"));
    }

    public void onError(Throwable err) throws Throwable {
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
            queryManager = new QueryManager(message, queryLogger.getLogger());
            if (queryManager.dryRun) {
                dryRunCycles = ((Double) queryManager.getQuery().func.aux.getOrDefault("cycles", 3.0)).intValue();
                logger.info(EventType.MISC.getMarker(), "Running dry run {} times", dryRunCycles);
                int dryRunResult = 0;
                int ordinaryResult = 0;
                for (int i = 0; i < dryRunCycles; i++) {
                    DataFrame dryRunDf = queryManager.run(true, false);
                    queryManager.closeConnection();
                    DataFrame ordinaryDf = queryManager.run(true, true);
                    queryManager.closeConnection();
                    if (i > 0) {
                        dryRunResult += (int) dryRunDf.columns.get(0).get(0);
                        ordinaryResult += (int) ordinaryDf.columns.get(0).get(0);
                    }
                }
                dryRunDiff = (ordinaryResult / dryRunCycles) - (dryRunResult / dryRunCycles);
                session.getRemote().sendStringByFuture("COMPLETED 0");
                return;
            }
            queryManager.initResultSet();
            if (queryManager.isResultSetInitialized()) {
                queryManager.initScheme();
                dataFrame = queryManager.getSubDF(dfNumber);
            } else {
                dataFrame = new DataFrame();
            }
        } else if (message.startsWith(SIZE_RECIEVED_MESSAGE)) {
            fdf = threadPool.submit(() -> queryManager.getSubDF(dfNumber + 1));
            Marker start = EventType.DATA_SENDING.getMarker(dfNumber, EventType.Stage.START);
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
            Marker start = EventType.DATAFRAME_TO_BYTEARRAY_CONVERTING.getMarker(dfNumber, EventType.Stage.START);
            Marker finish = EventType.DATAFRAME_TO_BYTEARRAY_CONVERTING.getMarker(dfNumber, EventType.Stage.END);
            logger.debug(start, "Converting dataframe to byteArray");
            bytes = dataFrame.toByteArray();
            logger.debug(finish, "DataFrame with id {} was converted to byteArray", dfNumber);
            logger.debug(EventType.CHECKSUM_SENDING.getMarker(dfNumber, EventType.Stage.START), "Sending checkSum message with bytes length of {}", bytes.length);
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
        logger.debug(EventType.LOG_SENDING.getMarker(EventType.Stage.START), "Converting logs to binary data and sending them to the server");
        DataFrame logs = queryLogger.dumpLogMessages();
        if (queryManager.dryRun) {
            logs.addRow(System.currentTimeMillis() * 1000.0, "GrokConnect", Level.INFO.levelStr,
                    EventType.DRY_RUN_DIFF.getName(), EventType.Stage.END.toString(), null, null, null, String.format("Average difference in duration when read result set and not for %s cycles", dryRunCycles - 1), dryRunDiff);
        }
        session.getRemote().sendBytesByFuture(ByteBuffer.wrap(logs.toByteArray()));
    }
}

package grok_connect.handlers;

import com.google.gson.Gson;
import grok_connect.connectors_info.DataQueryRunResult;
import grok_connect.log.EventType;
import grok_connect.log.QueryLogger;
import grok_connect.utils.QueryChunkNotSent;
import grok_connect.utils.QueryManager;
import org.eclipse.jetty.websocket.api.Session;
import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.Queue;
import java.util.concurrent.*;
import java.util.stream.Collectors;
import grok_connect.GrokConnect;
import org.slf4j.Logger;
import org.slf4j.Marker;
import serialization.DataFrame;

public class SessionHandler {
    private static final String COMPLETED = "COMPLETED_OK";
    private static final String MESSAGE_START = "QUERY";
    private static final String OK_RESPONSE = "DATAFRAME PART OK";
    private static final String END_MESSAGE = "EOF";
    private static final String SIZE_RECIEVED_MESSAGE = "DATAFRAME PART SIZE RECEIVED";
    private final Queue<Integer> chunkNumberQueue;
    private final Session session;
    private final ExecutorService threadPool;
    private final QueryLogger<DataFrame> queryLogger;
    private final Logger logger;
    private Future<DataFrame> fdf;
    private DataFrame dataFrame;
    private Boolean firstTry = true;
    private Boolean oneDfSent = false;
    private int messageNumber = 1;
    private byte[] bytes;
    private QueryManager queryManager;

    SessionHandler(Session session, QueryLogger<DataFrame> queryLogger) {
        threadPool = Executors.newCachedThreadPool();
        chunkNumberQueue = new LinkedBlockingQueue<>();
        this.queryLogger = queryLogger;
        this.logger = queryLogger.getLogger();
        this.session = session;
        logger.info(EventType.MISC.getMarker(), "GrokConnect version: {}", GrokConnect.properties.getProperty("version"));
    }

    public void onError(Throwable err) throws Throwable {
        Throwable cause = err.getCause() == null ? err : err.getCause(); // we need to get cause, because error is wrapped by runtime exception
        if (cause instanceof OutOfMemoryError) {
            // guess it won't work because there is no memory left!
            logger.error(EventType.ERROR.getMarker(), "Out of memory", cause);
            GrokConnect.needToReboot = true;
        } else {
            String message = cause.getMessage();
            String stackTrace = Arrays.stream(cause.getStackTrace()).map(StackTraceElement::toString)
                    .collect(Collectors.joining(System.lineSeparator()));
            logger.error(EventType.ERROR.getMarker(), "An exception was thrown", cause);
            DataQueryRunResult result = new DataQueryRunResult();
            result.errorMessage = message;
            result.errorStackTrace = stackTrace;
            result.log = queryLogger.dumpLogMessages().toCsv();
            session.getRemote().sendString(socketErrorMessage(new Gson().toJson(result)));
            session.close();
            queryManager.closeConnection();
        }
    }

    public void onMessage(String message) throws Throwable {
        logger.debug(EventType.SOCKET_MESSAGE_PROCESSING.getMarkerNumbered(messageNumber), "Received message {}", message);
        if (message.startsWith(MESSAGE_START)) {
            message = message.substring(6);
            queryManager = new QueryManager(message, queryLogger.getLogger());
            queryManager.initResultSet();
            if (queryManager.isResultSetInitialized()) {
                queryManager.initScheme();
                dataFrame = queryManager.getSubDF();
            } else {
                dataFrame = new DataFrame();
            }
            chunkNumberQueue.add(Integer.valueOf(dataFrame.tags.get(QueryManager.CHUNK_NUMBER_TAG)));
        } else if (message.startsWith(SIZE_RECIEVED_MESSAGE)) {
            long startTime = System.currentTimeMillis();
            fdf = threadPool.submit(() -> queryManager.getSubDF());
            Integer dfNumber = chunkNumberQueue.poll();
            Marker markerNumbered = EventType.CHUNK_SENDING.getMarkerNumbered(dfNumber);
            logger.debug(markerNumbered, "Sending binary dataframe with id {}", dfNumber);
            session.getRemote().sendBytes(ByteBuffer.wrap(bytes));
            logger.debug(markerNumbered, "Binary dataframe with id {} was sent to the server, exec time: {} ms", dfNumber,
                    System.currentTimeMillis() - startTime);
            return;
        } else if (message.startsWith(COMPLETED)) {
            logger.debug(EventType.LOG_PROCESSING.getMarker(), "Converting logs to byteArray");
            byte[] logs = queryLogger.dumpLogMessages().toByteArray();
            logger.debug(EventType.LOG_PROCESSING.getMarker(), "Logs were converted to byteArray");
            logger.debug(EventType.LOG_PROCESSING.getMarker(), "Sending logs");
            session.getRemote().sendBytes(ByteBuffer.wrap(logs));
            queryLogger.closeLogger();
            session.getRemote().sendString(END_MESSAGE);
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
                    chunkNumberQueue.add(Integer.valueOf(dataFrame.tags.get(QueryManager.CHUNK_NUMBER_TAG)));
                }
            }
        }
        if (dataFrame != null && (dataFrame.rowCount != 0 || !oneDfSent)) {
            Integer peek = chunkNumberQueue.peek();
            Marker marker = EventType.DATAFRAME_TO_BYTEARRAY_CONVERTING.getMarkerNumbered(peek);
            logger.debug(marker, "Converting dataframe to byteArray");
            long start = System.currentTimeMillis();
            bytes = dataFrame.toByteArray();
            logger.debug(marker, "DataFrame with id {} was converted to byteArray, execution time: {} ms", peek,
                    System.currentTimeMillis() - start);
            logger.debug(EventType.CHECKSUM_SENDING.getMarkerNumbered(peek), "Sending checkSum message with bytes length of {}", bytes.length);
            session.getRemote().sendString(checksumMessage(bytes.length));
        } else {
            logger.debug(EventType.MISC.getMarker(), "Closing connection");
            queryManager.closeConnection();
            logger.debug(EventType.MISC.getMarker(), "DB connection was closed");
            session.getRemote().sendString("COMPLETED");
        }
        logger.debug(EventType.SOCKET_MESSAGE_PROCESSING.getMarkerNumbered(messageNumber++), "Message was proceeded");
    }

    public String checksumMessage(int i) {
        return String.format("DATAFRAME PART SIZE: %d", i);
    }

    public QueryManager getQueryManager() {
        return queryManager;
    }

    private String socketErrorMessage(String s) {
        return String.format("ERROR: %s", s);
    }
}

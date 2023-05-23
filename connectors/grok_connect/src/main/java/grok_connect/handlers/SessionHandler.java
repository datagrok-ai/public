package grok_connect.handlers;

import com.google.gson.Gson;
import grok_connect.connectors_info.DataQueryRunResult;
import grok_connect.log.QueryLogger;
import grok_connect.utils.QueryChunkNotSent;
import grok_connect.utils.QueryManager;
import org.eclipse.jetty.websocket.api.Session;
import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.Collectors;
import grok_connect.GrokConnect;
import org.slf4j.Logger;
import serialization.DataFrame;

public class SessionHandler {
    private static final int rowsPerChunk = 10000;
    private static final String COMPLETED = "COMPLETED_OK";
    private static final String MESSAGE_START = "QUERY";
    private static final String OK_RESPONSE = "DATAFRAME PART OK";
    private static final String END_MESSAGE = "EOF";
    private static final String SIZE_RECIEVED_MESSAGE = "DATAFRAME PART SIZE RECEIVED";
    private final Session session;
    private final ExecutorService threadPool;
    private Future<DataFrame> fdf;
    private DataFrame dataFrame;
    private Boolean firstTry = true;
    private Boolean oneDfSent = false;
    private byte[] bytes;
    private QueryManager queryManager;
    private QueryLogger queryLogger;
    private Logger logger;

    SessionHandler(Session session) {
        threadPool = Executors.newCachedThreadPool();
        this.session = session;
    }

    public void onError(Throwable err) throws Throwable {
        Throwable cause = err.getCause(); // we need to get cause, because error is wrapped by runtime exception
        if (cause instanceof OutOfMemoryError) {
            // guess it won't work because there is no memory left!
            logger.error("Out of memory");
            GrokConnect.needToReboot = true;
        } else {
            String message = cause.getMessage();
            String stackTrace = Arrays.stream(cause.getStackTrace()).map(StackTraceElement::toString)
                    .collect(Collectors.joining(System.lineSeparator()));
            logger.warn("An exception was thrown");
            logger.debug("Error message: {}", message);
            logger.debug("Stack trace: {}", stackTrace);
            DataQueryRunResult result = new DataQueryRunResult();
            result.errorMessage = message;
            result.errorStackTrace = stackTrace;
//            result.log = queryLogger.dumpLogMessages();
            session.getRemote().sendString(socketErrorMessage(new Gson().toJson(result)));
            session.close();
            queryManager.closeConnection();
        }

    }

    public void onMessage(String message) throws Throwable {
        if (message.startsWith(MESSAGE_START)) {
            message = message.substring(6);
            queryManager = new QueryManager(message);
            queryLogger = queryManager.getQueryLogger();
            logger = queryLogger.getLogger();
            logger.info("GrokConnect version: {}", GrokConnect.properties.getProperty("version"));
            logger.debug("Received message that starts with {}", MESSAGE_START);
            queryManager.initResultSet();
            if (queryManager.isResultSetInitialized()) {
                queryManager.initScheme();
                dataFrame = queryManager.getSubDF(100);
            } else {
                dataFrame = new DataFrame();
            }
        } else if (message.startsWith(SIZE_RECIEVED_MESSAGE)) {
            logger.debug("Received message that starts with {}", SIZE_RECIEVED_MESSAGE);
            fdf = threadPool.submit(() -> queryManager.getSubDF(rowsPerChunk));
            logger.debug("Sending bytes");
            session.getRemote().sendBytes(ByteBuffer.wrap(bytes));
            return;
        } else if (message.startsWith(COMPLETED)) {
            logger.debug("Received message that starts with {}", COMPLETED);
            session.getRemote().sendBytes(ByteBuffer.wrap(queryLogger.dumpLogMessages().toByteArray()));
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
                logger.debug("Received message that starts with {}", OK_RESPONSE);
                firstTry = true;
                oneDfSent = true;
                if (queryManager.isResultSetInitialized())
                    dataFrame = fdf.get();
            }
        }
        if (dataFrame != null && (dataFrame.rowCount != 0 || !oneDfSent)) {
            long start = System.currentTimeMillis();
            bytes = dataFrame.toByteArray();
            logger.debug("DataFrame was converted to byteArray, execution time: {} ms",
                    System.currentTimeMillis() - start);
            logger.debug("Sending checkSum message with bytes length of {}", bytes.length);
            session.getRemote().sendString(checksumMessage(bytes.length));
        } else {
            logger.debug("Closing connection");
            queryManager.closeConnection();
            logger.debug("DB connection was closed");
            session.getRemote().sendString("COMPLETED");
        }
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

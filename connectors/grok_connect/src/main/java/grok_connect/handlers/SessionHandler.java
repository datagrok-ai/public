package grok_connect.handlers;

import com.google.gson.reflect.TypeToken;
import grok_connect.log.QueryLogger;
import grok_connect.utils.QueryChunkNotSent;
import grok_connect.utils.QueryManager;
import org.eclipse.jetty.websocket.api.Session;
import com.google.gson.Gson;
import java.nio.ByteBuffer;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import grok_connect.GrokConnect;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.DataFrame;

public class SessionHandler {
    private static final Logger LOGGER = LoggerFactory.getLogger(SessionHandler.class);
    private static final Gson gson = new Gson();
    private static final int rowsPerChunk = 10000;
    private static final String MESSAGE_START = "QUERY";
    private static final String OK_RESPONSE = "DATAFRAME PART OK";
    private static final String END_MESSAGE = "EOF";
    private static final String SIZE_RECIEVED_MESSAGE = "DATAFRAME PART SIZE RECEIVED";
    private static final String LOG_RECIEVED_MESSAGE = "LOG RECEIVED";
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
        dataFrame = null; // free some memory, maybe gc is benevolent today
        Throwable cause = err.getCause(); // we need to get cause, because error is wrapped by runtime exception
        String message = socketErrorMessage(cause);
        if (cause instanceof OutOfMemoryError) {
            // guess it won't work because there is no memory left!
            LOGGER.error(message);
            GrokConnect.needToReboot = true;
        } else {
            LOGGER.debug(message);
        }
        session.getRemote().sendString(message);
        session.close();
        queryManager.closeConnection();
    }

    public void onMessage(String message) throws Throwable {
        if (message.startsWith(MESSAGE_START)) {
            message = message.substring(6);
            queryManager = new QueryManager(message);
            queryLogger = queryManager.getQueryLogger();
            logger = queryLogger.getLogger();
            logger.info("GrokConnect version: {}", GrokConnect.properties.getProperty("version"));
            logger.debug("Received messages that starts with {}", MESSAGE_START);
            queryManager.initResultSet();
            if (queryManager.isResultSetInitialized()) {
                queryManager.initScheme();
                dataFrame = queryManager.getSubDF(100);
            } else {
                dataFrame = new DataFrame();
            }

        } else if (message.startsWith(LOG_RECIEVED_MESSAGE)) {
            logger.debug("Received messages that starts with {}", LOG_RECIEVED_MESSAGE);
            session.getRemote().sendString(END_MESSAGE);
            session.close();
            return;
        } else if (message.startsWith(SIZE_RECIEVED_MESSAGE)) {
            logger.debug("Received messages that starts with {}", SIZE_RECIEVED_MESSAGE);
            fdf = threadPool.submit(() -> queryManager.getSubDF(rowsPerChunk));
            logger.debug("Sending bytes");
            session.getRemote().sendBytes(ByteBuffer.wrap(bytes));
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
                if (queryManager.isResultSetInitialized())
                    dataFrame = fdf.get();
            }
        }
        if (dataFrame != null && (dataFrame.rowCount != 0 || !oneDfSent)) {
            long start = System.currentTimeMillis();
            bytes = dataFrame.toByteArray();
            logger.debug("DataFrame was converted to byteArray, execution time: {} ms",
                    System.currentTimeMillis() - start);
            logger.debug("Sending checkSum message");
            session.getRemote().sendString(checksumMessage(bytes.length));
        } else {
            logger.debug("Closing connection");
            queryManager.closeConnection();
            session.getRemote().sendString(socketLogMessage(queryLogger.dumpLogMessages()));
            queryLogger.closeLogger();
        }
    }

    public String checksumMessage(int i) {
        return String.format("DATAFRAME PART SIZE: %d", i);
    }

    public String socketLogMessage(String s) {
        return String.format("LOG: %s", s);
    }

    public QueryManager getQueryManager() {
        return queryManager;
    }

    private String socketErrorMessage(Throwable th) {
        Map<String, String> stackTrace = GrokConnect.printError(th);
        stackTrace.put("log", queryLogger.dumpLogMessages());
        return String.format("ERROR: %s", gson.toJson(stackTrace,
                new TypeToken<Map<String, String>>() { }.getType()));
    }
}

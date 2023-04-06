package grok_connect.handlers;

import grok_connect.connectors_info.FuncCall;
import grok_connect.utils.QueryCancelledByUser;
import grok_connect.utils.QueryChunkNotSent;
import grok_connect.utils.QueryManager;
import org.eclipse.jetty.websocket.api.Session;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.reflect.TypeToken;
import java.nio.ByteBuffer;
import java.text.SimpleDateFormat;
import java.util.Date;
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
    private static final int rowsPerChunk = 10000;
    private static final String OK_RESPONSE = "DATAFRAME PART OK";
    private static final String END_MESSAGE = "EOF";
    private static final String SIZE_RECIEVED_MESSAGE = "DATAFRAME PART SIZE RECEIVED";
    private static final String LOG_RECIEVED_MESSAGE = "LOG RECEIVED";
    private static final String DEFAULT_GETTING_RESULT_SET_MESSAGE = "Getting resultSet, %s";
    private static final String DEFAULT_GOT_RESULT_SET_MESSAGE = "Got resultSet, %s";
    private static final String DEFAULT_INITIALIZING_SCHEME_MESSAGE = "Initializing scheme, %s";
    private static final String DEFAULT_FINISHED_INITIALIZATION_SCHEME_MESSAGE = "Finished initialization of scheme, %s";
    private static final SimpleDateFormat DEFAULT_DATE_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.SSS");
    private final Session session;
    private final ExecutorService threadPool;
    private Future<DataFrame> fdf;
    private DataFrame dataFrame;
    private Boolean firstTry = true;
    private Boolean oneDfSent = false;
    private byte[] bytes;
    private QueryManager queryManager;
    private int queryMessages = 0;

    SessionHandler(Session session) {
        threadPool = Executors.newCachedThreadPool();
        this.session = session;
    }

    public void onError(Throwable err) throws Throwable {
        session.getRemote().sendString(socketErrorMessage(err));
        session.close();
        if (queryManager != null)
            queryManager.closeConnection();
    }

    public void onMessage(String message) throws Throwable {
        try {
            if (message.startsWith("QUERY")) {
                message = message.substring(6);
                queryManager = new QueryManager(message);
                FuncCall query = queryManager.getQuery();
                if (query.debugQuery) {
                    query.log += getOnMessageLogString(DEFAULT_GETTING_RESULT_SET_MESSAGE);
                }
                queryManager.initResultSet();
                if (query.debugQuery) {
                    query.log += getOnMessageLogString(DEFAULT_GOT_RESULT_SET_MESSAGE);
                }
                if (queryManager.isResultSetInitialized()) {
                    if (query.debugQuery) {
                        query.log += getOnMessageLogString(DEFAULT_INITIALIZING_SCHEME_MESSAGE);
                    }
                    queryManager.initScheme();
                    if (query.debugQuery) {
                        query.log += getOnMessageLogString(DEFAULT_FINISHED_INITIALIZATION_SCHEME_MESSAGE);
                    }
                    dataFrame = queryManager.getSubDF(100);
                } else {
                    dataFrame = new DataFrame();
                }

            } else if (message.startsWith(LOG_RECIEVED_MESSAGE)) {
                session.getRemote().sendString(END_MESSAGE);
                session.close();
                return;
            } else if (message.startsWith(SIZE_RECIEVED_MESSAGE)) {
                fdf = threadPool.submit(() -> queryManager.getSubDF(rowsPerChunk));
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
                bytes = dataFrame.toByteArray();
                session.getRemote().sendString(checksumMessage(bytes.length));
            } else {
                queryManager.closeConnection();
                FuncCall query = queryManager.getQuery();
                if (query.debugQuery) {
                    session.getRemote().sendString(socketLogMessage(query.log));
                    return;
                }

                session.getRemote().sendString(END_MESSAGE);
                session.close();
            }

        } catch (Throwable ex) {
            if (ex instanceof OutOfMemoryError) {
                LOGGER.error("Memory issue", ex);
                GrokConnect.needToReboot = true;
            }
            if (ex.getCause() instanceof QueryCancelledByUser) {
                LOGGER.debug("An exception was thrown", ex.getCause());
                ex = ex.getCause();
            }
            session.getRemote().sendString(socketErrorMessage(ex));
            session.close();
            queryManager.closeConnection();
        }

    }

    public String checksumMessage(int i) {
        return String.format("DATAFRAME PART SIZE: %d", i);
    }

    public String socketErrorMessage(Throwable th) {
        Gson gson = new GsonBuilder().create();
        return "ERROR: ".concat(gson.toJson(GrokConnect.printError(th),
                new TypeToken<Map<String, String>>() { }.getType()));
    }

    public String socketLogMessage(String s) {
        return "LOG: %s" + s;
    }

    public QueryManager getQueryManager() {
        return queryManager;
    }

    private String getOnMessageLogString(String messageFormat) {
        String log = String.format(messageFormat, DEFAULT_DATE_FORMAT.format(new Date()));
        LOGGER.trace(log);
        return log;
    }
}

package grok_connect.handlers;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import grok_connect.log.QueryLogger;
import grok_connect.log.QueryLoggerImpl;
import org.eclipse.jetty.websocket.api.Session;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class SessionManager {
    private static final String WRITE_LOG_HEADER = "writeLog";
    private static final Logger LOGGER = LoggerFactory.getLogger(SessionManager.class);
    private static final HashMap<Session, SessionHandler> sessions = new HashMap<>();
    private static final List<String> DEFAULT_LOG_LEVELS = Arrays.asList("info", "warn", "error");

    static void add(Session session) throws IOException {
        LOGGER.trace("add method was called with parameter: {}", session);
        String writeLog = session.getUpgradeRequest().getHeader(WRITE_LOG_HEADER);
        List<String> allowedLevels = new ArrayList<>();
        if (writeLog == null || writeLog.equals("false"))
            allowedLevels.addAll(DEFAULT_LOG_LEVELS);
        QueryLogger logger = new QueryLoggerImpl(session, allowedLevels);
        sessions.put(session, new SessionHandler(session, logger));
        logger.getLogger().debug("Sending CONNECTED message to the server...");
        session.getRemote().sendString("CONNECTED");
    }

    static void onMessage(Session session, String message) throws Throwable {
        LOGGER.trace("onMessage method was called with parameters: session :{}, message: {}", session, message);
        sessions.get(session).onMessage(message);
    }

    static void onError(Session session, Throwable error) throws Throwable {
        LOGGER.trace("onError method was called with parameters: session :{}, error: {}", session, error);
        sessions.get(session).onError(error);
    }

    static void delete(Session session) throws Throwable {
        LOGGER.trace("delete method was called with parameters: session :{}", session);
        if (sessions.containsKey(session)) {
            SessionHandler handler = sessions.get(session);
            handler.onClose();
            sessions.remove(session);
        }
    }
}

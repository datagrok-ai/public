package grok_connect.handlers;

import java.io.IOException;
import java.util.HashMap;
import grok_connect.utils.QueryManager;
import org.eclipse.jetty.websocket.api.Session;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class SessionManager {
    private static final Logger LOGGER = LoggerFactory.getLogger(SessionManager.class);
    private static final HashMap<Session, SessionHandler> sessions = new HashMap<>();

    static void add(Session session) throws IOException {
        LOGGER.trace("add method was called with parameter: {}", session);
        sessions.put(session, new SessionHandler(session));
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
            QueryManager queryManager = sessions.get(session).getQueryManager();
            if (queryManager != null)
                queryManager.closeConnection();
            sessions.remove(session);
        }
    }
}

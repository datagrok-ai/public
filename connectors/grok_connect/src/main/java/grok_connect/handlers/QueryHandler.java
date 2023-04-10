package grok_connect.handlers;

import java.io.IOException;
import org.eclipse.jetty.websocket.api.Session;
import org.eclipse.jetty.websocket.api.annotations.OnWebSocketClose;
import org.eclipse.jetty.websocket.api.annotations.OnWebSocketConnect;
import org.eclipse.jetty.websocket.api.annotations.OnWebSocketError;
import org.eclipse.jetty.websocket.api.annotations.OnWebSocketMessage;
import org.eclipse.jetty.websocket.api.annotations.WebSocket;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

@WebSocket
public class QueryHandler {
    private static final Logger LOGGER = LoggerFactory.getLogger(QueryHandler.class);

    @OnWebSocketConnect
    public void connected(Session session) throws IOException {
        LOGGER.trace("connected method was called with session: {}", session);
        SessionManager.add(session);
    }

    @OnWebSocketClose
    public void closed(Session session, int statusCode, String reason) throws Throwable {
        LOGGER.trace("closed method was called for session {} with statusCode {} and reason {}",
                session, statusCode, reason);
        SessionManager.delete(session);
    }

    @OnWebSocketMessage
    public void message(Session session, String message) throws Throwable {
        LOGGER.trace("recieved message {} for session {}", session, message);
        SessionManager.onMessage(session, message);
    }

    @OnWebSocketError
    public void error(Session session, Throwable error) throws Throwable {
        LOGGER.trace("an exception was thrown for {}. Exception: {}", session, error);
        SessionManager.onError(session, error);
    }
}

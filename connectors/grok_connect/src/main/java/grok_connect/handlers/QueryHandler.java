package grok_connect.handlers;

import java.io.IOException;
import org.eclipse.jetty.websocket.api.Session;
import org.eclipse.jetty.websocket.api.annotations.OnWebSocketClose;
import org.eclipse.jetty.websocket.api.annotations.OnWebSocketConnect;
import org.eclipse.jetty.websocket.api.annotations.OnWebSocketError;
import org.eclipse.jetty.websocket.api.annotations.OnWebSocketMessage;
import org.eclipse.jetty.websocket.api.annotations.WebSocket;

@WebSocket(maxTextMessageSize = 1024 * 1024) // set textMessageSize to 1MB
public class QueryHandler {

    @OnWebSocketConnect
    public void connected(Session session) throws IOException {
        session.setIdleTimeout(0);
        SessionManager.add(session);
    }

    @OnWebSocketClose
    public void closed(Session session, int statusCode, String reason) throws Throwable {
        SessionManager.delete(session);
    }

    @OnWebSocketMessage
    public void message(Session session, String message) throws Throwable {
        SessionManager.onMessage(session, message);
    }

    @OnWebSocketError
    public void error(Session session, Throwable error) {
        SessionManager.onError(session, error);
    }
}

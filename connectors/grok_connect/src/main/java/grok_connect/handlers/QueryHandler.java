package grok_connect.handlers;

import java.io.IOException;

import org.eclipse.jetty.websocket.api.*;
import org.eclipse.jetty.websocket.api.annotations.*;

// TODO add logger

@WebSocket
public class QueryHandler {

    @OnWebSocketConnect
    public void connected(Session session) throws IOException {
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
    public void error(Session session, Throwable error) throws Throwable {
        SessionManager.onError(session, error);
    }

}

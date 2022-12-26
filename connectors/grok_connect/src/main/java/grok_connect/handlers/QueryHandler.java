package grok_connect.handlers;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.sql.ResultSet;

import org.eclipse.jetty.websocket.api.*;
import org.eclipse.jetty.websocket.api.annotations.*;
import org.joda.time.*;
import com.google.gson.*;
import com.mysql.cj.xdevapi.Result;

import org.apache.log4j.*;

import serialization.*;
import grok_connect.utils.*;
import grok_connect.GrokConnect;
import grok_connect.connectors_info.*;
import grok_connect.providers.JdbcDataProvider;

// TODO add logger

@WebSocket
public class QueryHandler {

    QueryType qt;

    public QueryHandler(QueryType qt) {
        this.qt = qt;
    }

    @OnWebSocketConnect
    public void connected(Session session) throws IOException {
        SessionManager.add(session, qt);
    }

    @OnWebSocketClose
    public void closed(Session session, int statusCode, String reason) {
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

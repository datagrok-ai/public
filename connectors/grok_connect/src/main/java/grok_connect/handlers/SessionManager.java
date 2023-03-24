package grok_connect.handlers;

import java.io.IOException;
import java.util.HashMap;

import org.eclipse.jetty.websocket.api.Session;

public class SessionManager {
    static HashMap<Session, SessionHandler> sessions = new HashMap<Session, SessionHandler>();

    static void add(Session s) throws IOException {
        
        sessions.put(s, new SessionHandler(s));
        s.getRemote().sendString("CONNECTED");
    }

    static void onMessage(Session s, String message) throws Throwable {
        sessions.get(s).onMessage(message);
    }

    static void onError(Session s, Throwable error) throws Throwable {
        sessions.get(s).onError(error);
    }

    static void delete(Session s) throws Throwable {
        if (sessions.containsKey(s)) {
            SessionHandler sh = sessions.get(s);
            if (sh.qm != null)
                sessions.get(s).qm.closeConnection();
            sessions.remove(s);
        } 
    }
}

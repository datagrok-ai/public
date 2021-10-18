package grok_connect.utils;

import java.sql.Connection;
import java.sql.SQLException;
import java.util.*;

public class ConnectionPool {
    Map<String, Connection> connectionPull = Collections.synchronizedMap(new LinkedHashMap<>());
    int maxConnections = 100;

    ConnectionPool() {
        TimerTask timerTask = new TimerTask() {
            @Override
            public void run() {
                if (connectionPull.size() > maxConnections) {
                    for (Object conn : Arrays.asList(connectionPull.keySet().toArray()).subList(0, connectionPull.size() - maxConnections)) {
                        try {
                            connectionPull.get(conn).close();
                        } catch (SQLException throwables) {
                            throwables.printStackTrace(System.out);
                        }
                        connectionPull.remove(conn);
                    }
                }
            }
        };
        Timer timer = new Timer("MyTimer");
        timer.scheduleAtFixedRate(timerTask, 1*6*1000, 1*6*1000);
    }

    Connection getConnection(String url, java.util.Properties info, String driverClassName) {
        if (url != null && info != null && driverClassName != null && connectionPull.containsKey(url + info + driverClassName)) {
            Connection conn = connectionPull.get(url + info + driverClassName);
            try {
                if (!conn.isClosed())
                    return connectionPull.get(url + info + driverClassName);
            } catch (SQLException throwables) {
                //TODO: log in query
                throwables.printStackTrace(System.out);
            }
        }
        return null;
    }

    void putConnection(Connection conn, String url, java.util.Properties info, String driverClassName) {
        if (url != null && info != null && driverClassName != null)
            connectionPull.put(url + info + driverClassName, conn);
    }
}

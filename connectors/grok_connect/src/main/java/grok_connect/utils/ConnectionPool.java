package grok_connect.utils;

import com.zaxxer.hikari.HikariConfig;

import java.sql.Connection;
import java.sql.SQLException;
import java.util.*;

public class ConnectionPool {
    private static volatile ConnectionPool instance;
    private Timer timer;
    TimerTask timerTask;

    public static ConnectionPool getInstance() {
        ConnectionPool result = instance;
        if (result != null) {
            return result;
        }
        synchronized(ConnectionPool.class) {
            if (instance == null) {
                instance = new ConnectionPool();
            }
            return instance;
        }
    }

    public void setTimer() {

        if (timer == null)
            timer = new Timer("timer");

        if (timerTask != null)
            timerTask.cancel();

        timerTask = new TimerTask() {
            @Override
            public void run() {
                if (SettingsManager.getInstance().settings.debug) {
                    if (connectionPool.size() > 0)
                        System.out.println("Connection pool log:");
                    connectionPool.forEach((k, v) -> {
                        System.out.println("Pool: " + v.hikariDataSource.getPoolName());
                        System.out.println("Active: " + v.poolProxy.getActiveConnections());
                        System.out.println("Idle: " +v.poolProxy.getIdleConnections());
                    });
                }
            }
        };

        timer.scheduleAtFixedRate(timerTask,
                SettingsManager.getInstance().settings.connectionPoolTimerRate,
                SettingsManager.getInstance().settings.connectionPoolTimerRate);
    }

    public Map<String, HikariDataSourceInformation> connectionPool = Collections.synchronizedMap(new HashMap<>());

    public Connection getConnection(String url, java.util.Properties properties, String driverClassName) throws GrokConnectException, SQLException {
        if (url == null || properties == null || driverClassName == null)
            throw new GrokConnectException("Connection parameters are null");

        String key = url + properties + driverClassName;
        if (!connectionPool.containsKey(key))
            connectionPool.put(key, new HikariDataSourceInformation(url, properties, driverClassName));
        return connectionPool.get(key).hikariDataSource.getConnection();
    }
}

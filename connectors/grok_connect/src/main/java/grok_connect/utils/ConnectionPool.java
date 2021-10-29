package grok_connect.utils;

import com.zaxxer.hikari.HikariConfig;
import com.zaxxer.hikari.HikariDataSource;

import java.sql.Connection;
import java.sql.SQLException;
import java.util.*;

public class ConnectionPool {
    private static volatile ConnectionPool instance;

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

    public Map<String, HikariDataSource> connectionPool = Collections.synchronizedMap(new HashMap<>());

    public Connection getConnection(String url, java.util.Properties properties, String driverClassName) throws GrokConnectException, SQLException {
        if (url == null || properties == null || driverClassName == null)
            throw new GrokConnectException("Connection parameters are null");

        String key = url + properties + driverClassName;
        if (!connectionPool.containsKey(key)) {
            HikariConfig config = new HikariConfig();
            config.setJdbcUrl(url);
            config.setDataSourceProperties(properties);
            config.setDriverClassName(driverClassName);
            connectionPool.put(key, new HikariDataSource(config));
        }
        return connectionPool.get(key).getConnection();
    }
}

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

    Map<String, HikariDataSource> connectionPool = Collections.synchronizedMap(new LinkedHashMap<>());

    public Connection getConnection(String url, java.util.Properties properties, String driverClassName) {
        if (url != null && properties != null && driverClassName != null) {
            String key = url + properties + driverClassName;
            if (!connectionPool.containsKey(key)) {
                HikariConfig config = new HikariConfig();
                config.setJdbcUrl(url);
                config.setDataSourceProperties(properties);
                config.setDriverClassName(driverClassName);
                connectionPool.put(key, new HikariDataSource(config));
            }
            try {
                return connectionPool.get(key).getConnection();
            } catch (SQLException throwables) {
                //TODO: log in query
                throwables.printStackTrace(System.out);
            }
        }
        return null;
    }
}

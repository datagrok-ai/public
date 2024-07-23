package grok_connect.utils;

import com.zaxxer.hikari.HikariDataSource;
import com.zaxxer.hikari.pool.HikariPool;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.sql.Connection;
import java.sql.SQLException;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class ConnectionPool {
    private static final Logger LOGGER = LoggerFactory.getLogger(ConnectionPool.class);
    private static volatile ConnectionPool instance;
    private final Map<String, HikariDataSourceInformation> connectionPool;

    private ConnectionPool() {
        connectionPool = Collections.synchronizedMap(new HashMap<>());
    }

    public static synchronized ConnectionPool getInstance() {
        if (instance == null)
            instance = new ConnectionPool();
        return instance;
    }

    public Connection getConnection(String url, java.util.Properties properties, String driverClassName)
            throws GrokConnectException {
        try {
            LOGGER.debug("getConnection was called for driver {} with url {}", driverClassName, url);
            if (url == null || properties == null || driverClassName == null)
                throw new GrokConnectException("Connection parameters are null");
            String key = url + properties + driverClassName;
            synchronized(this) {
                if (!connectionPool.containsKey(key))
                    connectionPool.put(key, new HikariDataSourceInformation(url, properties, driverClassName));
                HikariDataSource hikariDataSource = connectionPool.get(key).hikariDataSource;
                return hikariDataSource.getConnection();
            }
        } catch (SQLException e) {
            throw new GrokConnectException(e);
        } catch (HikariPool.PoolInitializationException e) {
            throw new GrokConnectException(e.getCause());
        }
    }

    public Map<String, Connection> nativeConnectionsConnectionPool = Collections.synchronizedMap(new HashMap<>());

    Connection getNativeConnection(String url, java.util.Properties info, String driverClassName) {
        LOGGER.debug("getNativeConnection was called for driver {} with url {}", driverClassName, url);
        String key = url + info + driverClassName;
        if (url != null && info != null && driverClassName != null && nativeConnectionsConnectionPool.containsKey(key)) {
            Connection conn = nativeConnectionsConnectionPool.get(key);
            try {
                if (!conn.isClosed())
                    return nativeConnectionsConnectionPool.get(key);
            } catch (SQLException throwables) {
                LOGGER.warn("An exception was thrown", throwables);
            }
        }
        return null;
    }

    void putNativeConnection(Connection conn, String url, java.util.Properties info, String driverClassName) {
        if (url != null && info != null && driverClassName != null)
            nativeConnectionsConnectionPool.put(url + info + driverClassName, conn);
    }
}

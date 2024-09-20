package grok_connect.utils;

import com.zaxxer.hikari.pool.HikariPool;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.sql.Connection;
import java.sql.SQLException;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class ConnectionPool {
    private static final Logger LOGGER = LoggerFactory.getLogger(ConnectionPool.class);
    private static final Map<String, HikariDataSourceInformation> connectionPool = new ConcurrentHashMap<>();

    public static Connection getConnection(String url, java.util.Properties properties, String driverClassName)
            throws GrokConnectException {
        try {
            LOGGER.debug("getConnection was called for driver {} with url {}", driverClassName, url);
            if (GrokConnectUtil.isEmpty(url) || properties == null || GrokConnectUtil.isEmpty(driverClassName))
                throw new GrokConnectException("Connection parameters are null");
            String key = url + properties + driverClassName;
            HikariDataSourceInformation hikariDataSourceInformation = connectionPool.computeIfAbsent(key,
                    k -> new HikariDataSourceInformation(url, properties, driverClassName));
            return hikariDataSourceInformation.hikariDataSource.getConnection();
        } catch (SQLException e) {
            throw new GrokConnectException(e);
        } catch (HikariPool.PoolInitializationException e) {
            throw new GrokConnectException(e.getCause());
        }
    }
}

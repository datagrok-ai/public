package grok_connect.utils;

import com.zaxxer.hikari.HikariConfig;
import com.zaxxer.hikari.HikariDataSource;
import com.zaxxer.hikari.pool.HikariPool;
import grok_connect.connectors_info.DbCredentials;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.sql.DataSource;
import java.sql.Connection;
import java.sql.SQLException;
import java.sql.SQLTransientConnectionException;
import java.util.Map;
import java.util.Properties;
import java.util.concurrent.ConcurrentHashMap;

public class ConnectionPool {
    private static final Logger LOGGER = LoggerFactory.getLogger(ConnectionPool.class);
    private static final Map<String, HikariDataSource> connectionPool = new ConcurrentHashMap<>();
    private static final Object driverLoadLock = new Object();

    public static Connection getConnection(String url, java.util.Properties properties, String driverClassName)
            throws GrokConnectException {
        LOGGER.debug("getConnection was called for driver {} with url {}", driverClassName, url);
        if (GrokConnectUtil.isEmpty(url) || properties == null || GrokConnectUtil.isEmpty(driverClassName))
            throw new GrokConnectException("Connection parameters are null");
        String key = url + properties + driverClassName;
        try {
            HikariDataSource ds = getOrCreateDataSource(key, url, properties, driverClassName);
            return ds.getConnection();
        } catch (HikariPool.PoolInitializationException | SQLTransientConnectionException e) {
            HikariDataSource pool = connectionPool.remove(key);
            if (pool != null)
                pool.close();
            Throwable cause = e.getCause();
            throw new GrokConnectException(cause != null ? cause : e);
        } catch (SQLException | RuntimeException e) {
            throw new GrokConnectException(e);
        }
    }

    public static Connection getConnection(DataSource dataSource, String key)
            throws GrokConnectException {
        LOGGER.debug("getConnection was called for DataSource");
        if (dataSource == null || GrokConnectUtil.isEmpty(key))
            throw new GrokConnectException("Connection parameters are null");
        try {
            HikariDataSource ds = getOrCreateDataSource(key, dataSource);
            return ds.getConnection();
        } catch (HikariPool.PoolInitializationException | SQLTransientConnectionException e) {
            HikariDataSource pool = connectionPool.remove(key);
            if (pool != null)
                pool.close();
            Throwable cause = e.getCause();
            throw new GrokConnectException(cause != null ? cause : e);
        } catch (SQLException e) {
            throw new GrokConnectException(e);
        }
    }

    private static HikariDataSource getOrCreateDataSource(String key, String url, Properties properties, String driverClassName) {
        HikariDataSource ds = connectionPool.get(key);
        if (ds != null)
            return ds;
        HikariDataSource newDs = getDataSource(url, properties, driverClassName);
        HikariDataSource existing = connectionPool.putIfAbsent(key, newDs);
        if (existing != null) {
            newDs.close();
            return existing;
        }
        return newDs;
    }

    private static HikariDataSource getOrCreateDataSource(String key, DataSource dataSource) {
        HikariDataSource ds = connectionPool.get(key);
        if (ds != null)
            return ds;
        HikariDataSource newDs = getDataSource(dataSource);
        HikariDataSource existing = connectionPool.putIfAbsent(key, newDs);
        if (existing != null) {
            newDs.close();
            return existing;
        }
        return newDs;
    }

    private static HikariDataSource getDataSource(DataSource dataSource) {
        LOGGER.info("Initializing pool from DataSource");

        String poolName = "CustomDS - " +  " - " + dataSource.hashCode();
        LOGGER.debug("Pool name: {}", poolName);

        HikariConfig config = new HikariConfig();
        config.setPoolName(poolName);
        config.setDataSource(dataSource);
        config.setMaximumPoolSize(SettingsManager.getInstance().getSettings().connectionPoolMaximumPoolSize);
        config.setMinimumIdle(0);
        config.setIdleTimeout(SettingsManager.getInstance().getSettings().connectionPoolIdleTimeout);
        config.setLeakDetectionThreshold(60 * 1000);
        config.setInitializationFailTimeout(-1);

        return new HikariDataSource(config);
    }

    private static HikariDataSource getDataSource(String url, Properties properties, String driverClassName) {
        LOGGER.info("Initializing pool for driver {} with url {}", driverClassName, url);
        Properties propertiesWithoutPass = new Properties();
        propertiesWithoutPass.putAll(properties);
        propertiesWithoutPass.remove(DbCredentials.LOGIN);
        propertiesWithoutPass.remove(DbCredentials.PASSWORD);
        propertiesWithoutPass.remove(DbCredentials.ACCESS_KEY);
        propertiesWithoutPass.remove(DbCredentials.SECRET_KEY);
        propertiesWithoutPass.remove(DbCredentials.PRIVATE_KEY);
        propertiesWithoutPass.remove(DbCredentials.ACCOUNT_LOCATOR);
        propertiesWithoutPass.remove(DbCredentials.UID);
        propertiesWithoutPass.remove(DbCredentials.PWD);
        propertiesWithoutPass.remove(DbCredentials.OAUTH2_SECRET);
        propertiesWithoutPass.remove(DbCredentials.OAUTH2_CLIENT_ID);

        String alphanumeric = "[^A-Za-z\\d./|=]";
        String poolName = "Host - " + url.replaceAll("[:=]", "|").replaceAll(alphanumeric, "") +
                " . Properties - " + propertiesWithoutPass.toString().replaceAll("[:=]", "|").replaceAll(alphanumeric, "") +
                " . Driver - " + driverClassName.replaceAll("[:=]", "|").replaceAll(alphanumeric, "");
        LOGGER.debug("Pool name: {}", poolName);
        HikariConfig config = new HikariConfig();
        config.setPoolName(poolName);
        config.setJdbcUrl(url);
        synchronized (driverLoadLock) {
            config.setDriverClassName(driverClassName);
        }
        config.setDataSourceProperties(properties);
        config.setMaximumPoolSize(SettingsManager.getInstance().getSettings().connectionPoolMaximumPoolSize);
        config.setMinimumIdle(0);
        config.setIdleTimeout(SettingsManager.getInstance().getSettings().connectionPoolIdleTimeout);
        config.setLeakDetectionThreshold(60 * 1000);
        config.setInitializationFailTimeout(-1);
        return new HikariDataSource(config);
    }
}

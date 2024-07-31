package grok_connect.utils;

import com.zaxxer.hikari.HikariConfig;
import com.zaxxer.hikari.HikariDataSource;
import grok_connect.connectors_info.DbCredentials;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.util.Properties;

public class HikariDataSourceInformation {
    private static final Logger LOGGER = LoggerFactory.getLogger(HikariDataSourceInformation.class);
    public HikariDataSource hikariDataSource;

    public HikariDataSourceInformation(String url, Properties properties, String driverClassName) {
        LOGGER.info("Initializing pool for driver {} with url {}", driverClassName, url);
        Properties propertiesWithoutPass = new Properties();
        propertiesWithoutPass.putAll(properties);
        propertiesWithoutPass.remove(DbCredentials.LOGIN);
        propertiesWithoutPass.remove(DbCredentials.PASSWORD);
        propertiesWithoutPass.remove(DbCredentials.ACCESS_KEY);
        propertiesWithoutPass.remove(DbCredentials.SECRET_KEY);
        propertiesWithoutPass.remove(DbCredentials.ACCOUNT_LOCATOR);

        String alphanumeric = "[^A-Za-z0-9./|=]";
        String poolName = "Host - " + url.replaceAll("[:=]", "|").replaceAll(alphanumeric, "") +
                " . Properties - " + propertiesWithoutPass.toString().replaceAll("[:=]", "|").replaceAll(alphanumeric, "") +
                " . Driver - " + driverClassName.replaceAll("[:=]", "|").replaceAll(alphanumeric, "");
        LOGGER.debug("Pool name: {}", poolName);
        HikariConfig config = new HikariConfig();
        config.setPoolName(poolName);
        config.setJdbcUrl(url);
        config.setDriverClassName(driverClassName);
        config.setDataSourceProperties(properties);
        config.setMaximumPoolSize(SettingsManager.getInstance().getSettings().connectionPoolMaximumPoolSize);
        config.setMinimumIdle(0);
        config.setIdleTimeout(SettingsManager.getInstance().getSettings().connectionPoolIdleTimeout);
        config.setRegisterMbeans(false);

        this.hikariDataSource = new HikariDataSource(config);
    }
}

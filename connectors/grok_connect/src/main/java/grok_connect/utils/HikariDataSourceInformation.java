package grok_connect.utils;

import com.zaxxer.hikari.HikariConfig;
import com.zaxxer.hikari.HikariDataSource;
import com.zaxxer.hikari.HikariPoolMXBean;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import javax.management.JMX;
import javax.management.MBeanServer;
import javax.management.MalformedObjectNameException;
import javax.management.ObjectName;
import java.lang.management.ManagementFactory;
import java.util.Properties;

public class HikariDataSourceInformation {
    private static final Logger LOGGER = LoggerFactory.getLogger(HikariDataSourceInformation.class);
    public HikariDataSource hikariDataSource;
    public HikariPoolMXBean poolProxy;

    public HikariDataSourceInformation(String url, Properties properties, String driverClassName) {
        LOGGER.info("Initializing pool for driver {} with url {}", driverClassName, url);
        Properties propertiesWithoutPass = new Properties();
        propertiesWithoutPass.putAll(properties);
        propertiesWithoutPass.remove("password");

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
        config.setRegisterMbeans(true);

        this.hikariDataSource = new HikariDataSource(config);

        MBeanServer mBeanServer = ManagementFactory.getPlatformMBeanServer();
        ObjectName fullPoolName;
        try {
            fullPoolName = new ObjectName(String.format("com.zaxxer.hikari:type=Pool (%s)", poolName));
            LOGGER.debug("Fool pool name: {}", fullPoolName);
        } catch (MalformedObjectNameException e) {
            LOGGER.warn("An exception was thrown when creating fool pool name", e);
            throw new RuntimeException(e);
        }
        this.poolProxy = JMX.newMXBeanProxy(mBeanServer, fullPoolName, HikariPoolMXBean.class);
    }
}

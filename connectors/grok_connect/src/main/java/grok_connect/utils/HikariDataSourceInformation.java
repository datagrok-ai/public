package grok_connect.utils;

import com.zaxxer.hikari.HikariConfig;
import com.zaxxer.hikari.HikariDataSource;
import com.zaxxer.hikari.HikariPoolMXBean;

import javax.management.JMX;
import javax.management.MBeanServer;
import javax.management.MalformedObjectNameException;
import javax.management.ObjectName;
import java.lang.management.ManagementFactory;
import java.util.Properties;

public class HikariDataSourceInformation {
    public HikariDataSource hikariDataSource;
    public HikariPoolMXBean poolProxy;

    public HikariDataSourceInformation(String url, Properties properties, String driverClassName) {

        Properties propertiesWithoutPass = new Properties();
        propertiesWithoutPass.putAll(properties);
        propertiesWithoutPass.remove("password");

        String alphanumeric = "[^A-Za-z0-9./|=]";
        String poolName = "Host - " + url.replaceAll("[:=]", "|").replaceAll(alphanumeric, "") +
                " . Properties - " + propertiesWithoutPass.toString().replaceAll("[:=]", "|").replaceAll(alphanumeric, "") +
                " . Driver - " + driverClassName.replaceAll("[:=]", "|").replaceAll(alphanumeric, "");

        HikariConfig config = new HikariConfig();
        config.setPoolName(poolName);
        config.setJdbcUrl(url);
        config.setDriverClassName(driverClassName);
        config.setDataSourceProperties(properties);
        config.setMaximumPoolSize(SettingsManager.getInstance().settings.connectionPoolMaximumPoolSize);
        config.setMinimumIdle(0);
        config.setIdleTimeout(SettingsManager.getInstance().settings.connectionPoolIdleTimeout);
        config.setRegisterMbeans(true);

        this.hikariDataSource = new HikariDataSource(config);

        MBeanServer mBeanServer = ManagementFactory.getPlatformMBeanServer();
        ObjectName fullPoolName = null;
        try {
            fullPoolName = new ObjectName("com.zaxxer.hikari:type=Pool ("+poolName+")");
        } catch (MalformedObjectNameException e) {
            e.printStackTrace();
        }

        this.poolProxy = JMX.newMXBeanProxy(mBeanServer, fullPoolName, HikariPoolMXBean.class);
    }
}

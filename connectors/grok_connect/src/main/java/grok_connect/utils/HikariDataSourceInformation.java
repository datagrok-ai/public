package grok_connect.utils;

import com.zaxxer.hikari.HikariConfig;
import com.zaxxer.hikari.HikariDataSource;
import com.zaxxer.hikari.HikariPoolMXBean;

import javax.management.JMX;
import javax.management.MBeanServer;
import javax.management.MalformedObjectNameException;
import javax.management.ObjectName;
import java.lang.management.ManagementFactory;

public class HikariDataSourceInformation {
    public HikariDataSource hikariDataSource;

    public HikariDataSourceInformation(HikariConfig config) {
        config.setMaximumPoolSize(50);
        config.setMinimumIdle(0);
        config.setIdleTimeout(5*60*1000);

        this.hikariDataSource = new HikariDataSource(config);
    }
}

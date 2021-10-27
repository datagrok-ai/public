package grok_connect.providers;

import java.sql.*;
import java.util.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class MariaDbDataProvider extends MySqlDataProvider {
    public MariaDbDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "org.mariadb.jdbc.Driver";

        descriptor.type = "MariaDB";
        descriptor.description = "Query MariaDB database";
    }

    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString()) {
            properties.setProperty("zeroDateTimeBehavior", "convertToNull");
            if (conn.ssl()) {
                properties.setProperty("useSSL", "true");
                properties.setProperty("trustServerCertificate", "true");
            }
        }
        return properties;
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:mariadb://" + conn.getServer() + port + "/" + conn.getDb();
    }
}

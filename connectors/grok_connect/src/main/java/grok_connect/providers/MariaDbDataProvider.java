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

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName(driverClassName);
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString()) {
            properties.setProperty("zeroDateTimeBehavior", "convertToNull");
            if (conn.ssl()) {
                properties.setProperty("useSSL", "true");
                properties.setProperty("trustServerCertificate", "true");
            }
        }
        return CustomDriverManager.getConnection(getConnectionString(conn), properties, driverClassName);
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:mariadb://" + conn.getServer() + port + "/" + conn.getDb();
    }
}

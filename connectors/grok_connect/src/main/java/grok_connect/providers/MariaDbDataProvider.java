package grok_connect.providers;

import java.sql.*;
import java.util.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class MariaDbDataProvider extends MySqlDataProvider {
    public MariaDbDataProvider() {
        super();
        descriptor.type = "MariaDB";
        descriptor.description = "Query MariaDB database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.canBrowseSchema = true;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("org.mariadb.jdbc.Driver");
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString()) {
            properties.setProperty("zeroDateTimeBehavior", "convertToNull");
            if (conn.ssl()) {
                properties.setProperty("useSSL", "true");
                properties.setProperty("trustServerCertificate", "true");
            }
        }
        return DriverManager.getConnection(getConnectionString(conn), properties);
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:mariadb://" + conn.getServer() + port + "/" + conn.getDb();
    }
}

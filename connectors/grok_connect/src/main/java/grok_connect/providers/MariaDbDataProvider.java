package grok_connect.providers;

import java.util.Properties;
import grok_connect.connectors_info.DataConnection;
import grok_connect.resultset.DefaultResultSetManager;
import grok_connect.resultset.ResultSetManager;

public class MariaDbDataProvider extends MySqlDataProvider {
    public MariaDbDataProvider() {
        driverClassName = "org.mariadb.jdbc.Driver";
        descriptor.type = "MariaDB";
        descriptor.description = "Query MariaDB database";
    }

    @Override
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

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:mariadb://" + conn.getServer() + port + "/" + conn.getDb();
    }

    @Override
    public ResultSetManager getResultSetManager() {
        return DefaultResultSetManager.getDefaultManager();
    }
}

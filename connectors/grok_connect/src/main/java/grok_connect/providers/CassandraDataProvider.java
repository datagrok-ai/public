package grok_connect.providers;

import java.sql.*;
import java.util.*;
import org.apache.log4j.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class CassandraDataProvider extends JdbcDataProvider {
    public CassandraDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "com.github.cassandra.jdbc.CassandraDriver";

        descriptor = new DataSource();
        descriptor.type = "Cassandra";
        descriptor.description = "Query Cassandra database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.aggregations = null;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName(driverClassName);

        Logger logger = Logger.getLogger("com.github.cassandra.jdbc");
        logger.setLevel(Level.ERROR);

        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString() && conn.ssl()) {
            properties.setProperty("SSLMode", "1");
            properties.setProperty("UseSslIdentityCheck", "0");
        }
        return CustomDriverManager.getConnection(getConnectionString(conn), properties, driverClassName);
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:c*://" + conn.getServer() + port + "/" + conn.getDb();
    }
}

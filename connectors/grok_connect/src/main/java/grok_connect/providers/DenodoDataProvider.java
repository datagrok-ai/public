package grok_connect.providers;

import java.sql.*;
import java.util.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class DenodoDataProvider extends JdbcDataProvider {
    public DenodoDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "com.denodo.vdb.jdbcdriver.VDBJDBCDriver";

        descriptor = new DataSource();
        descriptor.type = "Denodo";
        descriptor.description = "Query Denodo database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName(driverClassName);
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString() && conn.ssl())
            properties.setProperty("ssl", "true");
        return CustomDriverManager.getConnection(getConnectionString(conn), properties, driverClassName);
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:vdb://" + conn.getServer() + port +  "/" + conn.getDb();
    }
}

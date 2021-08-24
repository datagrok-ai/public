package grok_connect.providers;

import java.sql.*;
import java.util.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class Db2DataProvider extends JdbcDataProvider {
    public Db2DataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "com.ibm.db2.jcc.DB2Driver";

        descriptor = new DataSource();
        descriptor.type = "DB2";
        descriptor.description = "Query DB2 database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName(driverClassName);
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString() && conn.ssl())
            properties.setProperty("sslConnection", "true");
        return CustomDriverManager.getConnection(getConnectionString(conn), properties, driverClassName);
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:db2://" + conn.getServer() + port + "/" + conn.getDb();
    }
}

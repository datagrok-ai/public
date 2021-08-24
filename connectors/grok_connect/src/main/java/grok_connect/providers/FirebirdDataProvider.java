package grok_connect.providers;

import java.sql.*;
import grok_connect.connectors_info.*;
import grok_connect.utils.CustomDriverManager;
import grok_connect.utils.ProviderManager;


public class FirebirdDataProvider extends JdbcDataProvider {
    public FirebirdDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "org.firebirdsql.jdbc.FBDriver";

        descriptor = new DataSource();
        descriptor.type = "Firebird";
        descriptor.description = "Query Firebird database";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName(driverClassName);
        return CustomDriverManager.getConnection(getConnectionString(conn), conn.credentials.getLogin(), conn.credentials.getPassword(), driverClassName);
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : "/" + conn.getPort();
        return "jdbc:firebirdsql:" + conn.getServer() + port + ":" + conn.getDb();
    }
}

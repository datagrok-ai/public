package grok_connect.providers;

import grok_connect.connectors_info.*;


public class FirebirdDataProvider extends JdbcDataProvider {
    public FirebirdDataProvider() {
        driverClassName = "org.firebirdsql.jdbc.FBDriver";

        descriptor = new DataSource();
        descriptor.type = "Firebird";
        descriptor.description = "Query Firebird database";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : "/" + conn.getPort();
        return "jdbc:firebirdsql:" + conn.getServer() + port + ":" + conn.getDb();
    }
}

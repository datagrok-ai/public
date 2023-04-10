package grok_connect.providers;

import grok_connect.connectors_info.*;
import grok_connect.resultset.ResultSetManager;
import grok_connect.utils.ProviderManager;


public class FirebirdDataProvider extends JdbcDataProvider {
    public FirebirdDataProvider(ResultSetManager resultSetManager, ProviderManager providerManager) {
        super(resultSetManager, providerManager);
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

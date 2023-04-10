package grok_connect.providers;

import grok_connect.connectors_info.*;
import grok_connect.resultset.ResultSetManager;
import grok_connect.utils.ProviderManager;


// TODO Note that this is just template. This class can be used as connector layer for other databases.
public class OdbcDataProvider extends JdbcDataProvider {
    public OdbcDataProvider(ResultSetManager resultSetManager, ProviderManager providerManager) {
        super(resultSetManager, providerManager);
        driverClassName = "sun.jdbc.odbc.JdbcOdbcDriver";

        descriptor = new DataSource();
        descriptor.type = "ODBC";
        descriptor.description = "Query database via ODBC";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:odbc://" + conn.getServer() + port + "/" + conn.getDb();
    }
}

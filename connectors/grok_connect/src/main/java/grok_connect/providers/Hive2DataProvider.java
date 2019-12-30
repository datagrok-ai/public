package grok_connect.providers;

import grok_connect.connectors_info.*;


public class Hive2DataProvider extends HiveDataProvider {
    public Hive2DataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Hive2";
        descriptor.description = "Query Hive2 database";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public String getConnectionString(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:hive2://" + conn.getServer() + port + "/" + conn.getDb();
    }
}

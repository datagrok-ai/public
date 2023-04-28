package grok_connect.providers;

import grok_connect.connectors_info.DataConnection;
import grok_connect.utils.ProviderManager;

public class Hive2DataProvider extends HiveDataProvider {
    public Hive2DataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "org.apache.hive.jdbc.HiveDriver";
        descriptor.type = "Hive2";
        descriptor.description = "Query Hive2 database";
    }

    @Override
    protected Integer getTimeout() {
        return 0;
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:hive2://" + conn.getServer() + port + "/" + conn.getDb();
    }
}

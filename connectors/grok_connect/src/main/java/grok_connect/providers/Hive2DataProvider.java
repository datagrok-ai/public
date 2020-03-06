package grok_connect.providers;

import java.util.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class Hive2DataProvider extends HiveDataProvider {
    public Hive2DataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Hive2";
        descriptor.description = "Query Hive2 database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:hive2://" + conn.getServer() + port + "/" + conn.getDb();
    }
}

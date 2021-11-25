package grok_connect.providers;

import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.table_query.AggrFunctionInfo;
import grok_connect.table_query.Stats;
import grok_connect.utils.Property;
import grok_connect.utils.ProviderManager;
import serialization.Types;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Properties;


public class СlickHouseProvider extends JdbcDataProvider {
    public СlickHouseProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "ru.yandex.clickhouse.ClickHouseDriver";

        descriptor = new DataSource();
        descriptor.type = "СlickHouse";
        descriptor.description = "Query СlickHouse database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:clickhouse://" + conn.getServer() + port + "/" + conn.getDb();
    }
}

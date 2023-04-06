package grok_connect.utils;

import grok_connect.connectors_info.DataSource;
import grok_connect.providers.AccessDataProvider;
import grok_connect.providers.AthenaDataProvider;
import grok_connect.providers.BigQueryDataProvider;
import grok_connect.providers.CassandraDataProvider;
import grok_connect.providers.ClickHouseProvider;
import grok_connect.providers.Db2DataProvider;
import grok_connect.providers.DenodoDataProvider;
import grok_connect.providers.DynamoDBDataProvider;
import grok_connect.providers.FirebirdDataProvider;
import grok_connect.providers.HBaseDataProvider;
import grok_connect.providers.Hive2DataProvider;
import grok_connect.providers.HiveDataProvider;
import grok_connect.providers.ImpalaDataProvider;
import grok_connect.providers.JdbcDataProvider;
import grok_connect.providers.MariaDbDataProvider;
import grok_connect.providers.MongoDbDataProvider;
import grok_connect.providers.MsSqlDataProvider;
import grok_connect.providers.MySqlDataProvider;
import grok_connect.providers.Neo4jDataProvider;
import grok_connect.providers.NeptuneDataProvider;
import grok_connect.providers.OracleDataProvider;
import grok_connect.providers.PIDataProvider;
import grok_connect.providers.PostgresDataProvider;
import grok_connect.providers.RedshiftDataProvider;
import grok_connect.providers.SQLiteDataProvider;
import grok_connect.providers.SapHanaDataProvider;
import grok_connect.providers.SnowflakeDataProvider;
import grok_connect.providers.TeradataDataProvider;
import grok_connect.providers.VerticaDataProvider;
import grok_connect.providers.VirtuosoDataProvider;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

public class ProviderManager {
    private static final Logger LOGGER = LoggerFactory.getLogger(ProviderManager.class);
    private final QueryMonitor queryMonitor;
    private final Map<String, JdbcDataProvider> providersMap;

    public ProviderManager() {
        this.queryMonitor = new QueryMonitor();
        LOGGER.debug("Initializing providers HashMap");
        List<JdbcDataProvider> providersList = new ArrayList<JdbcDataProvider>() {{
            add(new AccessDataProvider(ProviderManager.this));
            add(new AthenaDataProvider(ProviderManager.this));
            add(new BigQueryDataProvider(ProviderManager.this));
            add(new CassandraDataProvider(ProviderManager.this));
            add(new Db2DataProvider(ProviderManager.this));
            add(new FirebirdDataProvider(ProviderManager.this));
            add(new HBaseDataProvider(ProviderManager.this));
            add(new HiveDataProvider(ProviderManager.this));
            add(new Hive2DataProvider(ProviderManager.this));
            add(new MariaDbDataProvider(ProviderManager.this));
            add(new MongoDbDataProvider(ProviderManager.this));
            add(new MsSqlDataProvider(ProviderManager.this));
            add(new MySqlDataProvider(ProviderManager.this));
            add(new Neo4jDataProvider(ProviderManager.this));
            add(new OracleDataProvider(ProviderManager.this));
            add(new PostgresDataProvider(ProviderManager.this));
            add(new PIDataProvider(ProviderManager.this));
            add(new RedshiftDataProvider(ProviderManager.this));
            add(new SQLiteDataProvider(ProviderManager.this));
            add(new TeradataDataProvider(ProviderManager.this));
            add(new VerticaDataProvider(ProviderManager.this));
            add(new VirtuosoDataProvider(ProviderManager.this));
            add(new ImpalaDataProvider(ProviderManager.this));
            add(new DenodoDataProvider(ProviderManager.this));
            add(new SnowflakeDataProvider(ProviderManager.this));
            add(new ClickHouseProvider(ProviderManager.this));
            add(new NeptuneDataProvider(ProviderManager.this));
            add(new DynamoDBDataProvider(ProviderManager.this));
            add(new SapHanaDataProvider(ProviderManager.this));
        }};
        providersMap = providersList.stream()
                .collect(Collectors.toMap(provider -> provider.descriptor.type,
                Function.identity()));
    }

    public Collection<String> getAllProvidersTypes() {
        LOGGER.trace("getAllProvidersTypes providers was called");
        return providersMap.keySet();
    }

    public JdbcDataProvider getByName(String name) {
        LOGGER.trace("getByName with argument {} was called", name);
        return providersMap.get(name);
    }

    public QueryMonitor getQueryMonitor() {
        LOGGER.trace("getQueryMonitor was called");
        return queryMonitor;
    }

    public Collection<DataSource> getAllDescriptors() {
        LOGGER.trace("getAllDescriptors was called");
        return providersMap.values().stream()
                .map(provider -> provider.descriptor)
                .collect(Collectors.toList());
    }
}

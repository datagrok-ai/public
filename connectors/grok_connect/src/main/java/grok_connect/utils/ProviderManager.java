package grok_connect.utils;

import grok_connect.connectors_info.DataSource;
import grok_connect.providers.*;
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
    private final Map<String, JdbcDataProvider> providersMap;

    public ProviderManager() {
        LOGGER.debug("Initializing providers HashMap");
        List<JdbcDataProvider> providersList = new ArrayList<JdbcDataProvider>() {{
            add(new AccessDataProvider());
            add(new AthenaDataProvider());
            add(new BigQueryDataProvider());
            add(new CassandraDataProvider());
            add(new Db2DataProvider());
            add(new FirebirdDataProvider());
            add(new HBaseDataProvider());
            add(new HiveDataProvider());
            add(new Hive2DataProvider());
            add(new MariaDbDataProvider());
            add(new MongoDbDataProvider());
            add(new MsSqlDataProvider());
            add(new MySqlDataProvider());
            add(new Neo4jDataProvider());
            add(new OracleDataProvider());
            add(new PostgresDataProvider());
            add(new PIDataProvider());
            add(new RedshiftDataProvider());
            add(new SQLiteDataProvider());
            add(new TeradataDataProvider());
            add(new VerticaDataProvider());
            add(new VirtuosoDataProvider());
            add(new ImpalaDataProvider());
            add(new DenodoDataProvider());
            add(new SnowflakeDataProvider());
            add(new ClickHouseProvider());
            add(new NeptuneDataProvider());
            add(new DynamoDBDataProvider());
            add(new SapHanaDataProvider());
            add(new DatabricksProvider());
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

    public Collection<DataSource> getAllDescriptors() {
        LOGGER.trace("getAllDescriptors was called");
        return providersMap.values().stream()
                .map(provider -> provider.descriptor)
                .collect(Collectors.toList());
    }
}

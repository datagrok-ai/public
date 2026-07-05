package grok_connect.utils;

import grok_connect.connectors_info.DataSource;
import grok_connect.providers.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

public class ProviderManager {
    private static final Logger LOGGER = LoggerFactory.getLogger(ProviderManager.class);
    /**
     * Read-oriented engines / federation layers that must never advertise write support even though
     * they are JDBC providers with auto-interpolation (connector-writes WO-4).
     */
    private static final Set<String> WRITE_DENYLIST = new HashSet<>(Arrays.asList(
            "Athena", "BigQuery", "Impala", "Hive", "Hive2", "Virtuoso", "Denodo", "Neptune", "PI"));
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
//            add(new SQLiteDataProvider());
            add(new TeradataDataProvider());
            add(new VerticaDataProvider());
            add(new VirtuosoDataProvider());
            add(new ImpalaDataProvider());
            add(new DenodoDataProvider());
            add(new SnowflakeDataProvider());
            add(new ClickHouseProvider());
            add(new NeptuneDataProvider());
//            add(new DynamoDBDataProvider());
            add(new SapHanaDataProvider());
            add(new DatabricksProvider());
        }};
        providersMap = providersList.stream()
                .collect(Collectors.toMap(provider -> provider.descriptor.type,
                Function.identity()));
        // Central supportsWrite default: every auto-interpolation JDBC provider not on the denylist can
        // execute prepared-statement mutations. Providers that set the flag explicitly keep their value.
        for (JdbcDataProvider provider : providersMap.values()) {
            DataSource descriptor = provider.descriptor;
            if (!descriptor.supportsWrite)
                descriptor.supportsWrite = provider.autoInterpolation() && !WRITE_DENYLIST.contains(descriptor.type);
            // bulk insert rides the same prepared-statement batching as writes (default loader);
            // providers with a native fast path (Postgres COPY) still advertise it here (connector-writes WO-5).
            if (!descriptor.supportsBulkInsert)
                descriptor.supportsBulkInsert = descriptor.supportsWrite;
        }
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

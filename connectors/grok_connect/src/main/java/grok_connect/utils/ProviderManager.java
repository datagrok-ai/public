package grok_connect.utils;

import grok_connect.connectors_info.DataProvider;
import grok_connect.providers.*;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;

public class ProviderManager {
    private final Logger logger; // we reuse this logger in providers, is this good practice?
    private final QueryMonitor queryMonitor;
    public List<DataProvider> providers;

    public ProviderManager(Logger logger) {
        this.queryMonitor = new QueryMonitor();
        this.logger = logger;

        providers = new ArrayList<DataProvider>() {{
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
            add(new SapHanaDataProvider(ProviderManager.this));
        }};
    }

    public List<String> getAllProvidersTypes() {
        List<String> types = new ArrayList<>();

        for (DataProvider provider : providers)
            types.add(provider.descriptor.type);

        return types;
    }
    public DataProvider getByName(String name) {
        DataProvider provider = providers.get(0);

        for (ListIterator<DataProvider> providers = this.providers.listIterator(); providers.hasNext(); ) {
            DataProvider tmp = providers.next();
            if (tmp.descriptor.type.equals(name)) {
                provider = tmp;
                break;
            }
        }

        return provider;
    }

    public QueryMonitor getQueryMonitor() {
        return queryMonitor;
    }

    public Logger getLogger() {
        return logger;
    }
}

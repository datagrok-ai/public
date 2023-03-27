package grok_connect.utils;

import grok_connect.connectors_info.DataProvider;
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
import org.apache.log4j.Logger;
import java.util.ArrayList;
import java.util.List;

public class ProviderManager {
    private final Logger logger; // we reuse this logger in providers, is this good practice?
    private final QueryMonitor queryMonitor;
    public List<JdbcDataProvider> providers;

    public ProviderManager(Logger logger) {
        this.queryMonitor = new QueryMonitor();
        this.logger = logger;

        providers = new ArrayList<JdbcDataProvider>() {{
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
    }

    public List<String> getAllProvidersTypes() {
        List<String> types = new ArrayList<>();

        for (DataProvider provider : providers)
            types.add(provider.descriptor.type);

        return types;
    }
    public JdbcDataProvider getByName(String name) {
        JdbcDataProvider provider = providers.get(0);

        for (JdbcDataProvider tmp : this.providers) {
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

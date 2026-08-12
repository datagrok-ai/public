package grok_connect.utils;

import grok_connect.connectors_info.DataSource;
import grok_connect.providers.JdbcDataProvider;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class ProviderManager {
    private static final Logger LOGGER = LoggerFactory.getLogger(ProviderManager.class);

    // Registration is reflective so an image shipped without a provider's driver jar
    // (or even without the provider class) degrades to a startup warning instead of
    // failing class linking; the driver probe below keeps such providers off /conn.
    private static final String[] PROVIDER_CLASSES = {
            "grok_connect.providers.AccessDataProvider",
            "grok_connect.providers.AthenaDataProvider",
            "grok_connect.providers.BigQueryDataProvider",
            "grok_connect.providers.CassandraDataProvider",
            "grok_connect.providers.Db2DataProvider",
            "grok_connect.providers.FirebirdDataProvider",
            "grok_connect.providers.HBaseDataProvider",
            "grok_connect.providers.HiveDataProvider",
            "grok_connect.providers.Hive2DataProvider",
            "grok_connect.providers.MariaDbDataProvider",
            "grok_connect.providers.MongoDbDataProvider",
            "grok_connect.providers.MsSqlDataProvider",
            "grok_connect.providers.MySqlDataProvider",
            "grok_connect.providers.Neo4jDataProvider",
            "grok_connect.providers.OracleDataProvider",
            "grok_connect.providers.PostgresDataProvider",
            "grok_connect.providers.PIDataProvider",
            "grok_connect.providers.RedshiftDataProvider",
            "grok_connect.providers.TeradataDataProvider",
            "grok_connect.providers.VerticaDataProvider",
            "grok_connect.providers.VirtuosoDataProvider",
            "grok_connect.providers.ImpalaDataProvider",
            "grok_connect.providers.DenodoDataProvider",
            "grok_connect.providers.SnowflakeDataProvider",
            "grok_connect.providers.ClickHouseProvider",
            "grok_connect.providers.NeptuneDataProvider",
            "grok_connect.providers.SapHanaDataProvider",
            "grok_connect.providers.DatabricksProvider",
    };

    private final Map<String, JdbcDataProvider> providersMap;

    public ProviderManager() {
        Set<String> allowed = parseAllowlist(System.getenv("GROK_CONNECT_PROVIDERS"));
        providersMap = new LinkedHashMap<>();
        for (String className : PROVIDER_CLASSES) {
            JdbcDataProvider provider;
            try {
                provider = (JdbcDataProvider) Class.forName(className).newInstance();
            }
            catch (Throwable e) {
                LOGGER.warn("Provider {} skipped: {}", className, e.toString());
                continue;
            }
            String type = provider.descriptor.type;
            if (allowed != null && !allowed.contains(type))
                continue;
            String driver = provider.getDriverClassName();
            if (driver != null && !isDriverPresent(driver)) {
                LOGGER.warn("Provider {} skipped: driver {} is not on the classpath", type, driver);
                continue;
            }
            providersMap.put(type, provider);
        }
        LOGGER.info("Registered providers: {}", providersMap.keySet());
    }

    private static Set<String> parseAllowlist(String value) {
        if (GrokConnectUtil.isEmpty(value) || value.trim().equals("*"))
            return null;
        return Arrays.stream(value.split(","))
                .map(String::trim)
                .filter(GrokConnectUtil::isNotEmpty)
                .collect(Collectors.toCollection(LinkedHashSet::new));
    }

    private static boolean isDriverPresent(String driverClassName) {
        try {
            Class.forName(driverClassName, false, ProviderManager.class.getClassLoader());
            return true;
        }
        catch (Throwable e) {
            return false;
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

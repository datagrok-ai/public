package grok_connect.utils;

import grok_connect.column.*;
import grok_connect.connectors_info.DataSource;
import grok_connect.converter.ConverterManager;
import grok_connect.converter.array.ArrayConverterManager;
import grok_connect.converter.bigint.BigIntConverterManager;
import grok_connect.converter.bitstring.BitStringConverterManager;
import grok_connect.converter.bool.BoolTypeConverterManager;
import grok_connect.converter.complex.ComplexTypeConverterManager;
import grok_connect.converter.float_type.FloatTypeConverterManager;
import grok_connect.converter.integer.IntegerTypeConverterManager;
import grok_connect.converter.string.StringTypeConverterManager;
import grok_connect.converter.time.TimeTypeConverterManager;
import grok_connect.converter.xml.XmlConverterManager;
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
import grok_connect.resultset.DefaultResultSetManager;
import grok_connect.resultset.ResultSetManager;
import grok_connect.type.TypeChecker;
import grok_connect.type.array.DefaultArrayTypeChecker;
import grok_connect.type.bigint.ClickHouseBigIntTypeChecker;
import grok_connect.type.bigint.DefaultBigIntTypeChecker;
import grok_connect.type.bigint.Neo4jBigIntTypeChecker;
import grok_connect.type.bigint.OracleSnowflakeBigIntTypeChecker;
import grok_connect.type.bitstring.DefaultBitStringTypeChecker;
import grok_connect.type.bool.DefaultBoolTypeChecker;
import grok_connect.type.bool.MySQLMsSqlBoolTypeChecker;
import grok_connect.type.complex.DefaultComplexTypeChecker;
import grok_connect.type.complex.Neo4jComplexTypeChecker;
import grok_connect.type.float_type.DefaultFloatTypeChecker;
import grok_connect.type.integer.DefaultIntegerTypeChecker;
import grok_connect.type.integer.Neo4jIntegerTypeChecker;
import grok_connect.type.integer.OracleSnowflakeIntegerTypeChecker;
import grok_connect.type.string.DefaultStringTypeChecker;
import grok_connect.type.time.DefaultTimeTypeChecker;
import grok_connect.type.xml.DefaultXmlTypeChecker;
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
    private final ResultSetManager defaultManager;

    public ProviderManager() {
        this.queryMonitor = new QueryMonitor();
        LOGGER.debug("Initializing providers HashMap");
        defaultManager = getDefaultManager();
        ResultSetManager neo4jManager = getNeo4jManager();
        ResultSetManager oracleSnowflakeManager = getOracleSnowflakeManager();
        ResultSetManager clickHouseManager = getClickHouseManager();
        ResultSetManager mySQLManager = getMySQLManager();
        List<JdbcDataProvider> providersList = new ArrayList<JdbcDataProvider>() {{
            add(new AccessDataProvider(defaultManager, ProviderManager.this));
            add(new AthenaDataProvider(defaultManager, ProviderManager.this));
            add(new BigQueryDataProvider(defaultManager, ProviderManager.this));
            add(new CassandraDataProvider(defaultManager, ProviderManager.this));
            add(new Db2DataProvider(defaultManager,ProviderManager.this));
            add(new FirebirdDataProvider(defaultManager, ProviderManager.this));
            add(new HBaseDataProvider(defaultManager, ProviderManager.this));
            add(new HiveDataProvider(defaultManager, ProviderManager.this));
            add(new Hive2DataProvider(defaultManager, ProviderManager.this));
            add(new MariaDbDataProvider(defaultManager, ProviderManager.this));
            add(new MongoDbDataProvider(neo4jManager, ProviderManager.this));
            add(new MsSqlDataProvider(mySQLManager,ProviderManager.this));
            add(new MySqlDataProvider(mySQLManager,ProviderManager.this));
            add(new Neo4jDataProvider(neo4jManager, ProviderManager.this));
            add(new OracleDataProvider(oracleSnowflakeManager, ProviderManager.this));
            add(new PostgresDataProvider(defaultManager, ProviderManager.this));
            add(new PIDataProvider(defaultManager, ProviderManager.this));
            add(new RedshiftDataProvider(defaultManager, ProviderManager.this));
            add(new SQLiteDataProvider(defaultManager, ProviderManager.this));
            add(new TeradataDataProvider(defaultManager, ProviderManager.this));
            add(new VerticaDataProvider(defaultManager, ProviderManager.this));
            add(new VirtuosoDataProvider(defaultManager, ProviderManager.this));
            add(new ImpalaDataProvider(defaultManager, ProviderManager.this));
            add(new DenodoDataProvider(defaultManager, ProviderManager.this));
            add(new SnowflakeDataProvider(oracleSnowflakeManager, ProviderManager.this));
            add(new ClickHouseProvider(clickHouseManager, ProviderManager.this));
            add(new NeptuneDataProvider(defaultManager, ProviderManager.this));
            add(new DynamoDBDataProvider(defaultManager, ProviderManager.this));
            add(new SapHanaDataProvider(defaultManager, ProviderManager.this));
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

    public ResultSetManager getDefaultManager() {
        if (defaultManager != null) return defaultManager;
        List<ConverterManager<?>> converterManagers = new ArrayList<>();
        converterManagers.add(new BigIntConverterManager(new DefaultBigIntTypeChecker()));
        converterManagers.add(new IntegerTypeConverterManager(new DefaultIntegerTypeChecker()));
        converterManagers.add(new ComplexTypeConverterManager(new DefaultComplexTypeChecker()));
        initDefaults(converterManagers);
        return new DefaultResultSetManager(converterManagers, getDefaultColumnProviders());
    }

    private ResultSetManager getClickHouseManager() {
        List<ConverterManager<?>> converterManagers = new ArrayList<>();
        converterManagers.add(new BigIntConverterManager(new ClickHouseBigIntTypeChecker()));
        converterManagers.add(new IntegerTypeConverterManager(new DefaultIntegerTypeChecker()));
        converterManagers.add(new ComplexTypeConverterManager(new DefaultComplexTypeChecker()));
        initDefaults(converterManagers);
        return new DefaultResultSetManager(converterManagers, getClickHouseColumnProviders());
    }

    private ResultSetManager getOracleSnowflakeManager() {
        List<ConverterManager<?>> converterManagers = new ArrayList<>();
        converterManagers.add(new BigIntConverterManager(new OracleSnowflakeBigIntTypeChecker()));
        converterManagers.add(new IntegerTypeConverterManager(new OracleSnowflakeIntegerTypeChecker()));
        converterManagers.add(new ComplexTypeConverterManager(new DefaultComplexTypeChecker()));
        initDefaults(converterManagers);
        return new DefaultResultSetManager(converterManagers, getSnowflakeOracleColumnProviders());
    }

    private ResultSetManager getNeo4jManager() {
        List<ConverterManager<?>> converterManagers = new ArrayList<>();
        converterManagers.add(new BigIntConverterManager(new Neo4jBigIntTypeChecker()));
        converterManagers.add(new IntegerTypeConverterManager(new Neo4jIntegerTypeChecker()));
        converterManagers.add(new ComplexTypeConverterManager(new Neo4jComplexTypeChecker()));
        initDefaults(converterManagers);
        return new DefaultResultSetManager(converterManagers, getNeo4jColumnProviders());
    }

    private ResultSetManager getMySQLManager() {
        List<ConverterManager<?>> converterManagers = new ArrayList<>();
        converterManagers.add(new BigIntConverterManager(new DefaultBigIntTypeChecker()));
        converterManagers.add(new IntegerTypeConverterManager(new DefaultIntegerTypeChecker()));
        converterManagers.add(new ComplexTypeConverterManager(new DefaultComplexTypeChecker()));
        initDefaults(converterManagers);
        converterManagers.remove(5);
        converterManagers.add(new BoolTypeConverterManager(new MySQLMsSqlBoolTypeChecker()));
        return new DefaultResultSetManager(converterManagers, getMySQLColumnProviders());
    }

    private void initDefaults(List<ConverterManager<?>> converterManagers) {
        converterManagers.add(new ArrayConverterManager(new DefaultArrayTypeChecker()));
        converterManagers.add(new BitStringConverterManager(new DefaultBitStringTypeChecker()));
        converterManagers.add(new BoolTypeConverterManager(new DefaultBoolTypeChecker()));
        converterManagers.add(new FloatTypeConverterManager(new DefaultFloatTypeChecker()));
        converterManagers.add(new StringTypeConverterManager(new DefaultStringTypeChecker()));
        converterManagers.add(new TimeTypeConverterManager(new DefaultTimeTypeChecker()));
        converterManagers.add(new XmlConverterManager(new DefaultXmlTypeChecker()));
    }

    private Collection<TypeChecker> getStringColumnTypeCheckers() {
        List<TypeChecker> typeCheckers = new ArrayList<>();
        typeCheckers.add(new DefaultArrayTypeChecker());
        typeCheckers.add(new DefaultStringTypeChecker());
        typeCheckers.add(new DefaultXmlTypeChecker());
        return typeCheckers;
    }

    private Collection<TypeChecker> getBigIntColumnTypeCheckers() {
        List<TypeChecker> typeCheckers = new ArrayList<>();
        typeCheckers.add(new DefaultBigIntTypeChecker());
        typeCheckers.add(new DefaultBitStringTypeChecker());
        return typeCheckers;
    }

    private Collection<TypeChecker> getIntColumnTypeCheckers() {
        List<TypeChecker> typeCheckers = new ArrayList<>();
        typeCheckers.add(new DefaultIntegerTypeChecker());
        return typeCheckers;
    }

    private Collection<TypeChecker> getFloatColumnTypeCheckers() {
        List<TypeChecker> typeCheckers = new ArrayList<>();
        typeCheckers.add(new DefaultFloatTypeChecker());
        return typeCheckers;
    }

    private Collection<TypeChecker> getDateTimeColumnTypeCheckers() {
        List<TypeChecker> typeCheckers = new ArrayList<>();
        typeCheckers.add(new DefaultTimeTypeChecker());
        return typeCheckers;
    }

    private Collection<TypeChecker> getComplexColumnTypeCheckers() {
        List<TypeChecker> typeCheckers = new ArrayList<>();
        typeCheckers.add(new DefaultComplexTypeChecker());
        return typeCheckers;
    }

    private Collection<TypeChecker> getBoolColumnTypeCheckers() {
        List<TypeChecker> typeCheckers = new ArrayList<>();
        typeCheckers.add(new DefaultBoolTypeChecker());
        return typeCheckers;
    }

    private List<ColumnProvider> getCommonColumnProviders() {
        List<ColumnProvider> columnProviders = new ArrayList<>();
        ColumnProvider stringColumnProvider = new StringColumnProvider(getStringColumnTypeCheckers());
        columnProviders.add(stringColumnProvider);
        ColumnProvider floatColumnProvider = new FloatColumnProvider(getFloatColumnTypeCheckers());
        columnProviders.add(floatColumnProvider);
        ColumnProvider dateTimeColumnProvider = new DateTimeColumnProvider(getDateTimeColumnTypeCheckers());
        columnProviders.add(dateTimeColumnProvider);
        ColumnProvider boolColumnProvider = new BoolColumnProvider(getBoolColumnTypeCheckers());
        columnProviders.add(boolColumnProvider);
        return columnProviders;
    }

    private List<ColumnProvider> getDefaultColumnProviders() {
        List<ColumnProvider> columnProviders = new ArrayList<>();
        ColumnProvider intColumnProvider = new IntColumnProvider(getIntColumnTypeCheckers());
        columnProviders.add(intColumnProvider);
        ColumnProvider complexColumnProvider = new ComplexTypeColumnProvider(getComplexColumnTypeCheckers());
        columnProviders.add(complexColumnProvider);
        ColumnProvider bigIntColumnProvider = new BigIntColumnProvider(getBigIntColumnTypeCheckers());
        columnProviders.add(bigIntColumnProvider);
        columnProviders.addAll(getCommonColumnProviders());
        return columnProviders;
    }

    private List<ColumnProvider> getClickHouseColumnProviders() {
        List<ColumnProvider> columnProviders = new ArrayList<>();
        ColumnProvider intColumnProvider = new IntColumnProvider(getIntColumnTypeCheckers());
        columnProviders.add(intColumnProvider);
        ColumnProvider complexColumnProvider = new ComplexTypeColumnProvider(getComplexColumnTypeCheckers());
        columnProviders.add(complexColumnProvider);
        ColumnProvider bigIntColumnProvider = new BigIntColumnProvider(getClickHouseBigIntColumnTypeCheckers());
        columnProviders.add(bigIntColumnProvider);
        columnProviders.addAll(getCommonColumnProviders());
        return columnProviders;
    }

    private Collection<TypeChecker> getClickHouseBigIntColumnTypeCheckers() {
        List<TypeChecker> typeCheckers = new ArrayList<>();
        typeCheckers.add(new ClickHouseBigIntTypeChecker());
        return typeCheckers;
    }

    private List<ColumnProvider> getSnowflakeOracleColumnProviders() {
        List<ColumnProvider> columnProviders = new ArrayList<>();
        ColumnProvider intColumnProvider = new IntColumnProvider(getSnowflakeOracleIntegerColumnTypeCheckers());
        columnProviders.add(intColumnProvider);
        ColumnProvider complexColumnProvider = new ComplexTypeColumnProvider(getComplexColumnTypeCheckers());
        columnProviders.add(complexColumnProvider);
        ColumnProvider bigIntColumnProvider = new BigIntColumnProvider(getSnowflakeOracleBigIntColumnTypeCheckers());
        columnProviders.add(bigIntColumnProvider);
        columnProviders.addAll(getCommonColumnProviders());
        return columnProviders;
    }

    private Collection<TypeChecker> getSnowflakeOracleBigIntColumnTypeCheckers() {
        List<TypeChecker> typeCheckers = new ArrayList<>();
        typeCheckers.add(new OracleSnowflakeBigIntTypeChecker());
        return typeCheckers;
    }

    private Collection<TypeChecker> getSnowflakeOracleIntegerColumnTypeCheckers() {
        List<TypeChecker> typeCheckers = new ArrayList<>();
        typeCheckers.add(new OracleSnowflakeIntegerTypeChecker());
        return typeCheckers;
    }

    private List<ColumnProvider> getNeo4jColumnProviders() {
        List<ColumnProvider> columnProviders = new ArrayList<>();
        ColumnProvider intColumnProvider = new IntColumnProvider(getNeo4jIntegerColumnTypeCheckers());
        columnProviders.add(intColumnProvider);
        ColumnProvider complexColumnProvider = new ComplexTypeColumnProvider(getNeo4jComplexColumnTypeCheckers());
        columnProviders.add(complexColumnProvider);
        ColumnProvider bigIntColumnProvider = new BigIntColumnProvider(getNeo4jBigIntColumnTypeCheckers());
        columnProviders.add(bigIntColumnProvider);
        columnProviders.addAll(getCommonColumnProviders());
        return columnProviders;
    }

    private Collection<TypeChecker> getNeo4jBigIntColumnTypeCheckers() {
        List<TypeChecker> typeCheckers = new ArrayList<>();
        typeCheckers.add(new Neo4jBigIntTypeChecker());
        return typeCheckers;
    }

    private Collection<TypeChecker> getNeo4jIntegerColumnTypeCheckers() {
        List<TypeChecker> typeCheckers = new ArrayList<>();
        typeCheckers.add(new Neo4jIntegerTypeChecker());
        return typeCheckers;
    }

    private Collection<TypeChecker> getNeo4jComplexColumnTypeCheckers() {
        List<TypeChecker> typeCheckers = new ArrayList<>();
        typeCheckers.add(new Neo4jComplexTypeChecker());
        return typeCheckers;
    }

    private Collection<TypeChecker> getMyQLBoolColumnTypeCheckers() {
        List<TypeChecker> typeCheckers = new ArrayList<>();
        typeCheckers.add(new MySQLMsSqlBoolTypeChecker());
        return typeCheckers;
    }

    private List<ColumnProvider> getMySQLColumnProviders() {
        List<ColumnProvider> columnProviders = new ArrayList<>();
        ColumnProvider intColumnProvider = new IntColumnProvider(getIntColumnTypeCheckers());
        columnProviders.add(intColumnProvider);
        ColumnProvider complexColumnProvider = new ComplexTypeColumnProvider(getComplexColumnTypeCheckers());
        columnProviders.add(complexColumnProvider);
        ColumnProvider bigIntColumnProvider = new BigIntColumnProvider(getBigIntColumnTypeCheckers());
        columnProviders.add(bigIntColumnProvider);
        List<ColumnProvider> commonColumnProviders = getCommonColumnProviders();
        commonColumnProviders.set(3, new BoolColumnProvider(getMyQLBoolColumnTypeCheckers()));
        columnProviders.addAll(commonColumnProviders);
        return columnProviders;
    }
}

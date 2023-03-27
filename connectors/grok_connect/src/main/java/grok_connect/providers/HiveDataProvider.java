package grok_connect.providers;

import java.io.IOException;
import java.sql.*;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Properties;
import java.util.stream.Collectors;
import grok_connect.connectors_info.*;
import grok_connect.providers.proxy.HiveMetaDataProviderProxyProvider;
import grok_connect.table_query.AggrFunctionInfo;
import grok_connect.table_query.GroupAggregation;
import grok_connect.table_query.Stats;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.Prop;
import grok_connect.utils.Property;
import grok_connect.utils.ProviderManager;
import grok_connect.utils.QueryCancelledByUser;
import serialization.DataFrame;
import serialization.Types;

public class HiveDataProvider extends JdbcDataProvider {
    private static final List<String> AVAILABLE_META_STORES =
            Collections.unmodifiableList(Arrays.asList("MySQL", "Postgres", "Oracle", "MS SQL"));
    private static final String DEFAULT_META_STORE_DB = "metastore";

    public HiveDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "org.apache.hadoop.hive.jdbc.HiveDriver";

        descriptor = new DataSource();
        descriptor.type = "Hive";
        descriptor.description = "Query Hive database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        Property metaStores = new Property(Property.STRING_TYPE, DbCredentials.META_STORE,
                "Select database engine hosting Hive Metastore to enable schema browsing support");
        metaStores.choices = AVAILABLE_META_STORES;
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.connectionTemplate.add(metaStores);
        descriptor.connectionTemplate.add(new Property(Property.STRING_TYPE, DbCredentials.META_STORE_SERVER,
                "Hostname or IP address of server on which Hive Metastore is available"));
        descriptor.connectionTemplate.add(new Property(Property.INT_TYPE, DbCredentials.META_STORE_PORT,
                "Port of Hive Metastore"));
        descriptor.connectionTemplate.add(new Property(Property.STRING_TYPE, DbCredentials.META_STORE_DB,
                "Database name of metastore. Defaults to \"metastore\""));
        descriptor.connectionTemplate.add(new Property(Property.STRING_TYPE, DbCredentials.META_STORE_LOGIN,
                "Hive Metastore login"));
        descriptor.connectionTemplate.add(new Property(Property.STRING_TYPE, DbCredentials.META_STORE_PASSWORD,
                "Hive Metastore password",
                new Prop("password")));
        descriptor.canBrowseSchema = true;
        descriptor.typesMap = new HashMap<String, String>() {{
            put("smallint", Types.INT);
            put("tinyint", Types.INT);
            put("int", Types.INT);
            put("bigint", Types.BIG_INT);
            put("boolean", Types.BOOL);
            put("double", Types.FLOAT);
            put("float", Types.FLOAT);
            put("#decimal.*", Types.FLOAT);
            put("date", Types.DATE_TIME);
            put("timestamp", Types.DATE_TIME);
            put("#char.*", Types.STRING);
            put("#varchar.*", Types.STRING);
            put("string", Types.STRING);
            put("#array.*", Types.OBJECT);
            put("#map.*", Types.OBJECT);
            put("#struct.*", Types.OBJECT);
        }};
        descriptor.aggregations.add(new AggrFunctionInfo(Stats.STDEV, "stddev(#)", Types.dataFrameNumericTypes));
    }

    @Override
    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString() && conn.ssl())
            properties.setProperty("ssl", "true");
        return properties;
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:hive://" + conn.getServer() + port + "/" + conn.getDb();
    }

    @Override
    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException, GrokConnectException {
        prepareProvider();
        return DriverManager.getConnection(getConnectionString(conn), getProperties(conn));
    }

    @Override
    public DataFrame getSchemas(DataConnection connection) throws ClassNotFoundException, SQLException,
            ParseException, IOException, QueryCancelledByUser, GrokConnectException {
        return getProxyMetaStoreProvider(connection).getSchemas(prepareMetaStoreConnection(connection));
    }

    @Override
    public DataFrame getSchema(DataConnection connection, String schema, String table) throws
            ClassNotFoundException, SQLException, ParseException, IOException, QueryCancelledByUser,
            GrokConnectException {
        return getProxyMetaStoreProvider(connection)
                .getSchema(prepareMetaStoreConnection(connection), schema, table);
    }

    @Override
    protected String getRegexQuery(String columnName, String regexExpression) {
        return String.format("%s RLIKE '%s'", columnName, regexExpression);
    }

    @Override
    protected void appendQueryParam(DataQuery dataQuery, String paramName, StringBuilder queryBuffer) {
        FuncParam param = dataQuery.getParam(paramName);
        if (param.propertyType.equals("list")) {
            @SuppressWarnings("unchecked")
            List<String> values = ((ArrayList<String>) param.value);
            queryBuffer.append(values.stream().map(value -> "?").collect(Collectors.joining(", ")));
        } else {
            queryBuffer.append("?");
        }
    }

    @Override
    protected int setArrayParamValue(PreparedStatement statement, int n, FuncParam param) throws SQLException {
        @SuppressWarnings (value="unchecked")
        ArrayList<Object> lst = (ArrayList<Object>)param.value;
        if (lst == null || lst.size() == 0) {
            statement.setObject(n, null);
            return 0;
        }
        for (int i = 0; i < lst.size(); i++) {
            statement.setObject(n + i, lst.get(i));
        }
        return lst.size() - 1;
    }

    private JdbcDataProvider getProxyMetaStoreProvider(DataConnection connection) {
        return new HiveMetaDataProviderProxyProvider()
                .getProxy(providerManager, getMetaDataProvider(connection).getClass());
    }

    private JdbcDataProvider getMetaDataProvider(DataConnection connection) {
        if (connection.parameters.containsKey(DbCredentials.META_STORE)) {
            return providerManager.getByName((String) connection.parameters.get(DbCredentials.META_STORE));
        }
        throw new UnsupportedOperationException("Hive Metastore information should be provided for schema browsing");
    }

    private DataConnection prepareMetaStoreConnection(DataConnection connection) {
        Credentials credentials = new Credentials();
        credentials.parameters.put(DbCredentials.LOGIN, connection.get(DbCredentials.META_STORE_LOGIN));
        credentials.parameters.put(DbCredentials.PASSWORD, connection.get(DbCredentials.META_STORE_PASSWORD));
        DataConnection metaStoreConnection = new DataConnection();
        metaStoreConnection.credentials = credentials;
        String metaStoreDb = (String) connection.parameters.get(DbCredentials.META_STORE_DB);
        metaStoreConnection.parameters.put(DbCredentials.DB, metaStoreDb == null ? DEFAULT_META_STORE_DB : metaStoreDb);
        metaStoreConnection.parameters.put(DbCredentials.SERVER, connection.get(DbCredentials.META_STORE_SERVER));
        metaStoreConnection.parameters.put(DbCredentials.PORT, connection.parameters.get(DbCredentials.META_STORE_PORT));
        return metaStoreConnection;
    }
}

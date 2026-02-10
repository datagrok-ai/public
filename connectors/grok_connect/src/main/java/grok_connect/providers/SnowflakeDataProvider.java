package grok_connect.providers;

import grok_connect.managers.ColumnManager;
import grok_connect.managers.bigint_column.SnowflakeBigIntColumnManager;
import grok_connect.managers.integer_column.OracleSnowflakeIntColumnManager;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataQuery;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncParam;
import grok_connect.resultset.DefaultResultSetManager;
import grok_connect.resultset.ResultSetManager;
import grok_connect.table_query.AggrFunctionInfo;
import grok_connect.table_query.Stats;
import grok_connect.utils.*;
import net.snowflake.client.jdbc.SnowflakeBasicDataSource;
import net.snowflake.client.pooling.SnowflakeConnectionPoolDataSource;
import serialization.Types;

import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.*;

public class SnowflakeDataProvider extends JdbcDataProvider {
    private static final boolean CAN_BROWSE_SCHEMA = true;
    private static final String DEFAULT_SCHEMA = "PUBLIC";
    private static final String URL_PREFIX = "jdbc:snowflake://";
    private static final String URL_SEPARATOR = ".";
    private static final String SERVER = "snowflakecomputing.com";
    private static final String DRIVER_CLASS_NAME = "net.snowflake.client.jdbc.SnowflakeDriver";
    private static final String TYPE = "Snowflake";
    private static final String DESCRIPTION = "Query Snowflake database";
    private static final List<String> AVAILABLE_CLOUDS =
            Collections.unmodifiableList(Arrays.asList("aws", "azure", "gcp", "privatelink"));
    private static final String RSA_METHOD = "Username/Private key";

    public SnowflakeDataProvider() {
        init();
    }

    @Override
    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = super.getProperties(conn);
        if (!conn.hasCustomConnectionString()) {
            setIfNotNull(properties, DbCredentials.DB, conn.getDb());
            setIfNotNull(properties, DbCredentials.WAREHOUSE, conn.get(DbCredentials.WAREHOUSE));
            setIfNotNull(properties, DbCredentials.ACCOUNT, buildAccount(conn));
            String schema = conn.get(DbCredentials.SCHEMA);
            properties.setProperty(DbCredentials.SCHEMA, schema == null ? DEFAULT_SCHEMA : schema);
            setIfNotNull(properties, DbCredentials.ROLE, conn.get(DbCredentials.ROLE));
        }
        return properties;
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        return new StringBuilder(URL_PREFIX)
                .append(buildAccount(conn))
                .append(URL_SEPARATOR)
                .append(SERVER)
                .toString();
    }

    @Override
    public Connection getConnection(DataConnection conn) throws SQLException, GrokConnectException {
        prepareProvider();
        if (conn.credentials == null)
            throw new GrokConnectException("Credentials can't be null");
        String method = (String) conn.credentials.parameters.get("#chosen-auth-method");
        if (GrokConnectUtil.isNotEmpty(method) && method.equals(RSA_METHOD)) {
            SnowflakeBasicDataSource ds = getDataSource(conn);
            return ConnectionPool.getConnection(ds, getDataSourceKey(conn, ds));
        }
        else
            return ConnectionPool.getConnection(getConnectionString(conn), getProperties(conn), driverClassName);
    }

    @Override
    public String getSchemasSql(String db) {
        return "SELECT DISTINCT table_schema FROM information_schema.columns;";
    }

    @Override
    public String getSchemaSql(String db, String schema, String table, boolean includeKeyInfo) {
        boolean isEmptyDb = GrokConnectUtil.isEmpty(db);
        boolean isEmptySchema = GrokConnectUtil.isEmpty(schema);
        boolean isEmptyTable = GrokConnectUtil.isEmpty(table);
        String whereClause = String.format(" WHERE%s%s%s",
                isEmptyDb ? "" : String.format(" LOWER(c.table_catalog) = LOWER('%s')", db),
                isEmptySchema ? "" : String.format("%s c.table_schema = '%s'", isEmptyDb ? "" : " AND", schema),
                isEmptyTable ? "" : String.format("%s c.table_name = '%s'", isEmptyDb && isEmptySchema ? "" : " AND", table));
        return String.format("SELECT c.table_catalog as table_catalog, c.table_schema as table_schema, c.table_name as table_name, c.column_name as column_name, "
                        + "c.data_type as data_type, "
                        + "case t.table_type when 'VIEW' then 1 else 0 end as is_view FROM information_schema.columns c "
                        + "JOIN information_schema.tables t ON t.table_name = c.table_name AND t.table_schema = c.table_schema AND LOWER(t.table_catalog) = LOWER(c.table_catalog)%s " +
                        "ORDER BY c.table_name, c.ordinal_position;",
                isEmptyDb && isEmptySchema && isEmptyTable ? "" : whereClause);
    }

    @Override
    protected String getRegexQuery(String columnName, String regexExpression) {
        return String.format("%s REGEXP '%s'", columnName, regexExpression);
    }

    @Override
    protected void appendQueryParam(DataQuery dataQuery, String paramName, StringBuilder queryBuffer) {
        FuncParam param = dataQuery.getParam(paramName);
        if (param.propertyType.equals("list")) {
            queryBuffer.append("SELECT TRIM(VALUE) FROM TABLE(SPLIT_TO_TABLE(?, ','))");
        } else {
            queryBuffer.append("?");
        }
    }

    @Override
    protected int setArrayParamValue(PreparedStatement statement, int n, FuncParam param) throws SQLException {
        @SuppressWarnings("unchecked")
        List<String> list = ((ArrayList<String>) param.value);
        String values = String.join(",", list);
        statement.setObject(n, values);
        return 0;
    }

    @Override
    public ResultSetManager getResultSetManager() {
        Map<String, ColumnManager<?>> defaultManagersMap = DefaultResultSetManager.getDefaultManagersMap();
        defaultManagersMap.put(Types.INT, new OracleSnowflakeIntColumnManager());
        defaultManagersMap.put(Types.BIG_INT, new SnowflakeBigIntColumnManager());
        return DefaultResultSetManager.fromManagersMap(defaultManagersMap);
    }

    private void init() {
        driverClassName = DRIVER_CLASS_NAME;
        descriptor = new DataSource();
        descriptor.type = TYPE;
        descriptor.description = DESCRIPTION;
        descriptor.canBrowseSchema = CAN_BROWSE_SCHEMA;
        descriptor.supportCatalogs = true;
        descriptor.defaultSchema = DEFAULT_SCHEMA;
        Property cloudProviders = new Property(Property.STRING_TYPE, DbCredentials.CLOUD);
        cloudProviders.choices = AVAILABLE_CLOUDS;
        descriptor.connectionTemplate = new ArrayList<Property>(){{
            add(new Property(Property.STRING_TYPE, DbCredentials.ACCOUNT_LOCATOR));
            add(new Property(Property.STRING_TYPE, DbCredentials.REGION_ID));
            add(cloudProviders);
            add(new Property(Property.STRING_TYPE, DbCredentials.DB, DbCredentials.DB_DESCRIPTION));
            add(new Property(Property.STRING_TYPE, DbCredentials.WAREHOUSE));
            add(new Property(Property.STRING_TYPE, DbCredentials.ROLE));
            add(new Property(Property.STRING_TYPE, DbCredentials.CONNECTION_STRING,
                    DbCredentials.CONNECTION_STRING_DESCRIPTION, new Prop("textarea")));
        }};
        descriptor.credentialsTemplate = DbCredentials.getDbCredentialsTemplate();
        descriptor.credentialsTemplate.add(new Property(Property.STRING_TYPE, DbCredentials.PRIVATE_KEY, null, RSA_METHOD, new Prop("rsa"), ".pem,.der"));
        descriptor.credentialsTemplate.add(new Property(Property.STRING_TYPE, "passPhrase", "Passphrase for decrypting private key.", RSA_METHOD, new Prop("password")));
        descriptor.credentialsTemplate.stream().filter(p -> p.name.equals(DbCredentials.LOGIN)).forEach(p -> p.category = "Username/Password," + RSA_METHOD);
        descriptor.nameBrackets = "\"";

        descriptor.typesMap = new HashMap<String, String>() {{
            put("number", Types.FLOAT);
            put("float", Types.FLOAT);
            put("text", Types.STRING);
            put("boolean", Types.BOOL);
            put("date", Types.DATE_TIME);
            put("time", Types.DATE_TIME);
            put("#timestamp.*", Types.DATE_TIME);
            put("array", Types.LIST);
            put("object", Types.OBJECT);
            put("variant", Types.OBJECT);
            put("geography", Types.OBJECT);
            put("binary", Types.BLOB);
        }};
        descriptor.aggregations.add(new AggrFunctionInfo(Stats.STDEV, "stddev(#)", Types.dataFrameNumericTypes));
        descriptor.jdbcPropertiesTemplate = new ArrayList<Property>() {{
            add(new Property(Property.INT_TYPE, "loginTimeout",
                    "Maximum time to wait for login in seconds", new Prop()));
            add(new Property(Property.INT_TYPE, "networkTimeout",
                    "Socket read/write timeout in milliseconds", new Prop()));
            add(new Property(Property.INT_TYPE, "queryTimeout",
                    "Query execution timeout in seconds", new Prop()));
            add(new Property(Property.BOOL_TYPE, "CLIENT_SESSION_KEEP_ALIVE",
                    "Keep the session alive until explicitly closed", new Prop(), "keepAlive"));


            add(new Property(Property.INT_TYPE, "CLIENT_PREFETCH_THREADS",
                    "Number of threads to prefetch result chunks", new Prop(), "prefetchThreads"));
            add(new Property(Property.BOOL_TYPE, "CLIENT_DISABLE_INCIDENTS",
                    "Disable sending incidents to Snowflake", new Prop(), "disableIncidents"));
        }};
    }

    private String buildAccount(DataConnection conn) {
        StringBuilder builder = new StringBuilder(conn.get(DbCredentials.ACCOUNT_LOCATOR));
        if (GrokConnectUtil.isNotEmpty(conn.get(DbCredentials.REGION_ID)))
            builder.append(URL_SEPARATOR).append(conn.get(DbCredentials.REGION_ID));
        if (GrokConnectUtil.isNotEmpty(conn.get(DbCredentials.CLOUD)))
            builder.append(URL_SEPARATOR).append(conn.get(DbCredentials.CLOUD));
        return builder.toString();
    }

    private SnowflakeBasicDataSource getDataSource(DataConnection conn) throws GrokConnectException {
        try {
            SnowflakeBasicDataSource sfDataSource = new SnowflakeConnectionPoolDataSource();
            sfDataSource.setPrivateKey(KeyLoader.get((String) conn.credentials.parameters.get(DbCredentials.PRIVATE_KEY), (String) conn.credentials.parameters.get("passPhrase")));
            sfDataSource.setUser(conn.credentials.getLogin());
            if (!conn.hasCustomConnectionString()) {
                if (conn.getDb() != null)
                    sfDataSource.setDatabaseName(conn.getDb());
                if (conn.get(DbCredentials.WAREHOUSE) != null)
                    sfDataSource.setWarehouse(conn.get(DbCredentials.WAREHOUSE));
                if (conn.get(DbCredentials.SCHEMA) != null)
                    sfDataSource.setSchema(conn.get(DbCredentials.SCHEMA));
                if (conn.get(DbCredentials.ROLE) != null)
                    sfDataSource.setRole(conn.get(DbCredentials.ROLE));
            }
            sfDataSource.setUrl(getConnectionString(conn));
            return sfDataSource;
        } catch (Exception e) {
            throw new GrokConnectException(e);
        }
    }

    private String getDataSourceKey(DataConnection conn, SnowflakeBasicDataSource ds) throws GrokConnectException {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            String privateKeyContent = (String) conn.credentials.parameters.get(DbCredentials.PRIVATE_KEY);
            String privateKeyHash = privateKeyContent != null
                    ? Base64.getEncoder().encodeToString(digest.digest(privateKeyContent.trim().getBytes(StandardCharsets.UTF_8)))
                    : "";

            String rawKey = String.join("|",
                    "user=" + GrokConnectUtil.nullSafe(conn.credentials.getLogin()),
                    "account=" + GrokConnectUtil.nullSafe(buildAccount(conn)),
                    "db=" + GrokConnectUtil.nullSafe(conn.getDb()),
                    "schema=" + GrokConnectUtil.nullSafe(conn.get(DbCredentials.SCHEMA)),
                    "warehouse=" + GrokConnectUtil.nullSafe(conn.get(DbCredentials.WAREHOUSE)),
                    "role=" + GrokConnectUtil.nullSafe(conn.get(DbCredentials.ROLE)),
                    "url=" + GrokConnectUtil.nullSafe(ds.getUrl()),
                    "key=" + privateKeyHash
            );

            return Base64.getEncoder().encodeToString(rawKey.getBytes(StandardCharsets.UTF_8));
        } catch (NoSuchAlgorithmException e) {
            throw new GrokConnectException(e);
        }
    }
}

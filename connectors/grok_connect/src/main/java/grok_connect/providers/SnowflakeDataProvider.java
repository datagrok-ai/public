package grok_connect.providers;

import grok_connect.connectors_info.*;
import grok_connect.table_query.AggrFunctionInfo;
import grok_connect.table_query.Stats;
import grok_connect.utils.Prop;
import grok_connect.utils.Property;
import grok_connect.utils.ProviderManager;
import serialization.Types;

import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Properties;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class SnowflakeDataProvider extends JdbcDataProvider{
    private static final boolean CAN_BROWSE_SCHEMA = true;
    private static final String DEFAULT_SCHEMA = "PUBLIC";
    private static final String URL_PREFIX = "jdbc:snowflake://";
    private static final String URL_SEPARATOR = ".";
    private static final String SERVER = "snowflakecomputing.com";
    private static final String DRIVER_CLASS_NAME = "net.snowflake.client.jdbc.SnowflakeDriver";
    private static final String TYPE = "Snowflake";
    private static final String DESCRIPTION = "Query Snowflake database";
    private static final List<String> AVAILABLE_CLOUDS =
            Collections.unmodifiableList(Arrays.asList("aws", "azure", "gcp"));

    public SnowflakeDataProvider(ProviderManager providerManager) {
        super(providerManager);
        init();
    }

    @Override
    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString()) {
            properties.put(DbCredentials.DB, conn.getDb());
            properties.put(DbCredentials.WAREHOUSE, conn.get(DbCredentials.WAREHOUSE));
            properties.put(DbCredentials.ACCOUNT, buildAccount(conn));
            String schema = conn.get(DbCredentials.SCHEMA);
            properties.put(DbCredentials.SCHEMA, schema == null ? DEFAULT_SCHEMA : schema);
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
    public String getSchemasSql(String db) {
        return "SELECT DISTINCT table_schema FROM information_schema.columns;";
    }

    @Override
    public String getSchemaSql(String db, String schema, String table) {
        schema = schema != null ? schema : DEFAULT_SCHEMA;
        String query = "SELECT table_schema, table_name, column_name, data_type "
                + "FROM information_schema.columns "
                + "WHERE table_schema = '%s' ORDER BY table_name";
        return String.format(query, schema);
    }

    @Override
    protected boolean isInteger(int type, String typeName, int precision, int scale) {
        return typeName.equals("NUMBER") && precision < 10 && scale == 0;
    }

    @Override
    protected boolean isBigInt(int type, String typeName, int precision, int scale) {
        return typeName.equals("NUMBER") && precision >= 10 && scale == 0;
    }

    @Override
    protected boolean isFloat(int type, String typeName, int precision, int scale) {
        return typeName.equals("DOUBLE") || type == java.sql.Types.DOUBLE;
    }

    @Override
    protected String getRegexQuery(String columnName, String regexExpression) {
        return String.format("%s REGEXP '%s'", columnName, regexExpression);
    }

    @Override
    protected void appendQueryParam(DataQuery dataQuery, String paramName, StringBuilder queryBuffer) {
        FuncParam param = dataQuery.getParam(paramName);
        if (param.propertyType.equals("list")) {
            queryBuffer.append("SELECT TRIM(VALUE) FROM TABLE(strtok_split_to_table('?', ','))");
        } else {
            queryBuffer.append("?");
        }
    }

    @Override
    protected int setArrayParamValue(PreparedStatement statement, int n, FuncParam param) throws SQLException {
        @SuppressWarnings("unchecked")
        String[] values = ((ArrayList<String>) param.value).toArray(new String[0]);
        statement.setString(n, String.join(",", values ));
        return 0;
    }

    private void init() {
        driverClassName = DRIVER_CLASS_NAME;
        descriptor = new DataSource();
        descriptor.type = TYPE;
        descriptor.description = DESCRIPTION;
        descriptor.canBrowseSchema = CAN_BROWSE_SCHEMA;
        descriptor.defaultSchema = DEFAULT_SCHEMA;
        Prop notNullProperty = new Prop();
        notNullProperty.nullable = false;
        Property cloudProviders = new Property(Property.STRING_TYPE, DbCredentials.CLOUD, notNullProperty);
        cloudProviders.choices = AVAILABLE_CLOUDS;
        descriptor.connectionTemplate = new ArrayList<Property>(){{
            add(new Property(Property.STRING_TYPE, DbCredentials.ACCOUNT_LOCATOR, notNullProperty));
            add(new Property(Property.STRING_TYPE, DbCredentials.REGION_ID, notNullProperty));
            add(cloudProviders);
            add(new Property(Property.STRING_TYPE, DbCredentials.DB, DbCredentials.DB_DESCRIPTION));
            add(new Property(Property.STRING_TYPE, DbCredentials.WAREHOUSE));
            add(new Property(Property.STRING_TYPE, DbCredentials.CONNECTION_STRING,
                    DbCredentials.CONNECTION_STRING_DESCRIPTION, new Prop("textarea")));
            add(new Property(Property.BOOL_TYPE, DbCredentials.CACHE_SCHEMA));
            add(new Property(Property.BOOL_TYPE, DbCredentials.CACHE_RESULTS));
            add(new Property(Property.STRING_TYPE, DbCredentials.CACHE_INVALIDATE_SCHEDULE));
        }};
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.nameBrackets = "\"";

        //TODO: .*
        descriptor.typesMap = new HashMap<String, String>() {{
            put("number", Types.FLOAT);
            put("decimal", Types.FLOAT);
            put("numeric", Types.FLOAT);
            put("int", Types.INT);
            put("integer", Types.INT);
            put("bigint", Types.BIG_INT);
            put("smallint", Types.INT);
            put("float", Types.FLOAT);
            put("float4", Types.FLOAT);
            put("float8", Types.FLOAT);
            put("double", Types.FLOAT);
            put("double precision", Types.FLOAT);
            put("real", Types.FLOAT);
            put("varchar", Types.STRING);
            put("char", Types.STRING);
            put("character", Types.STRING);
            put("string", Types.STRING);
            put("text", Types.STRING);
            put("boolean", Types.BOOL);
            put("date", Types.DATE_TIME);
            put("datetime", Types.DATE_TIME);
            put("time", Types.DATE_TIME);
            put("timestamp", Types.DATE_TIME);
        }};
        descriptor.aggregations.add(new AggrFunctionInfo(Stats.STDEV, "stddev(#)", Types.dataFrameNumericTypes));
    }

    private String buildAccount(DataConnection conn) {
        return new StringBuilder(conn.get(DbCredentials.ACCOUNT_LOCATOR))
                .append(URL_SEPARATOR)
                .append(conn.get(DbCredentials.REGION_ID))
                .append(URL_SEPARATOR)
                .append(conn.get(DbCredentials.CLOUD))
                .toString();
    }
}

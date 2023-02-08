package grok_connect.providers;

import com.sun.tools.javac.util.List;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.table_query.AggrFunctionInfo;
import grok_connect.table_query.Stats;
import grok_connect.utils.Prop;
import grok_connect.utils.Property;
import grok_connect.utils.ProviderManager;
import serialization.Types;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Properties;

public class SnowflakeDataProvider extends JdbcDataProvider{
    private static final String URL_PREFIX = "jdbc:snowflake://";
    private static final String URL_SEPARATOR = ".";
    private static final String SERVER = "snowflakecomputing.com";
    private static final String DRIVER_CLASS_NAME = "net.snowflake.client.jdbc.SnowflakeDriver";
    private static final String TYPE = "Snowflake";
    private static final String DESCRIPTION = "Query Snowflake database";
    private static final List<String> AVAILABLE_CLOUDS = List.of("aws", "azure", "gcp");

    public SnowflakeDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = DRIVER_CLASS_NAME;
        descriptor = new DataSource();
        descriptor.type = TYPE;
        descriptor.description = DESCRIPTION;
        Prop prop = new Prop();
        prop.nullable = false;
        Property cloudProviders = new Property(Property.STRING_TYPE, DbCredentials.CLOUD, prop);
        cloudProviders.choices = AVAILABLE_CLOUDS;
        descriptor.connectionTemplate = new ArrayList<Property>(){{
            add(new Property(Property.STRING_TYPE, DbCredentials.ACCOUNT_LOCATOR, prop));
            add(new Property(Property.STRING_TYPE, DbCredentials.REGION_ID, prop));
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

    @Override
    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString()) {
            properties.put(DbCredentials.DB, conn.getDb());
            properties.put(DbCredentials.WAREHOUSE, conn.get(DbCredentials.WAREHOUSE));
            properties.put(DbCredentials.ACCOUNT, buildAccount(conn));
        }
        return properties;
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        return URL_PREFIX + buildAccount(conn)
                + URL_SEPARATOR + SERVER;
    }

    private String buildAccount(DataConnection conn) {
        return conn.get(DbCredentials.ACCOUNT_LOCATOR)
                + URL_SEPARATOR + conn.get(DbCredentials.REGION_ID) + URL_SEPARATOR
                + conn.get(DbCredentials.CLOUD);
    }
}

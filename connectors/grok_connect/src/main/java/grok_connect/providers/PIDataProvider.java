package grok_connect.providers;

import java.util.*;

import grok_connect.resultset.ResultSetManager;
import serialization.Types;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class PIDataProvider extends JdbcDataProvider {
    public PIDataProvider(ResultSetManager resultSetManager, ProviderManager providerManager) {
        super(resultSetManager, providerManager);
        driverClassName = "com.osisoft.jdbc.Driver";

        descriptor = new DataSource();
        descriptor.type = "PI";
        descriptor.description = "Query PI database";
        descriptor.connectionTemplate = new ArrayList<Property>() {{
                                                add(new Property(Property.STRING_TYPE, DbCredentials.ACCESS_SERVER));
                                                add(new Property(Property.STRING_TYPE, DbCredentials.SERVER));
                                                add(new Property(Property.STRING_TYPE, DbCredentials.DB, DbCredentials.DB_DESCRIPTION));
                                                add(new Property(Property.STRING_TYPE, DbCredentials.CONNECTION_STRING,
                                                        DbCredentials.CONNECTION_STRING_DESCRIPTION, new Prop("textarea")));
                                                add(new Property(Property.BOOL_TYPE, DbCredentials.CACHE_SCHEMA));
                                                add(new Property(Property.BOOL_TYPE, DbCredentials.CACHE_RESULTS));
                                                add(new Property(Property.STRING_TYPE, DbCredentials.CACHE_INVALIDATE_SCHEDULE));
                                            }};
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.nameBrackets = "\"";

        descriptor.canBrowseSchema = true;
        descriptor.defaultSchema = "public";
        descriptor.typesMap = new HashMap<String, String>() {{
            put("smallint", Types.INT);
            put("int", Types.INT);
            put("bigint", Types.BIG_INT);
            put("real", Types.FLOAT);
            put("double precision", Types.FLOAT);
            put("numeric", Types.FLOAT);
            put("#character.*", Types.STRING);
            put("#varchar.*", Types.STRING);
            put("text", Types.STRING);
        }};
    }

    public String getConnectionStringImpl(DataConnection conn) {
        return "jdbc:pioledb://" + conn.parameters.get(DbCredentials.ACCESS_SERVER) + "/Data Source=" + conn.getServer() + ";Initial Catalog=" + conn.getDb() + ";Integrated Security=SSPI";
    }

    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = defaultConnectionProperties(conn);
        return properties;
    }
}

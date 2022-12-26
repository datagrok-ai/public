package grok_connect.providers;

import java.sql.*;
import java.util.*;
import serialization.Types;
import grok_connect.utils.*;
import grok_connect.table_query.*;
import grok_connect.connectors_info.*;


public class PIDataProvider extends JdbcDataProvider {
    public PIDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "com.osisoft.jdbc.Driver";

        descriptor = new DataSource();
        descriptor.type = "PI";
        descriptor.description = "Query PI database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
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

    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = defaultConnectionProperties(conn);
        return properties;
    }
}

package grok_connect.providers;

import java.util.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;

public class DenodoDataProvider extends JdbcDataProvider {
    public DenodoDataProvider() {
        driverClassName = "com.denodo.vdp.jdbc.Driver";

        descriptor = new DataSource();
        descriptor.type = "Denodo";
        descriptor.description = "Query Denodo database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
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
        return "jdbc:denodo://" + conn.getServer() + port +  "/" + conn.getDb();
    }
}

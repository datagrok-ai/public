package grok_connect.providers;

import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.utils.Property;
import java.util.ArrayList;

public class OrientDbJdbcProvider extends JdbcDataProvider {
    public OrientDbJdbcProvider() {
        driverClassName = "com.orientechnologies.orient.jdbc.OrientJdbcDriver";

        descriptor = new DataSource();
        descriptor.type = "OrientDb";
        descriptor.description = "Query OrientDb";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        return String.format("jdbc:orient:remote:%s/%s", conn.getServer(), conn.getDb());
    }
}

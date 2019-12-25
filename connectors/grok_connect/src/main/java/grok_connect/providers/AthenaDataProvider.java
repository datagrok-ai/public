package grok_connect.providers;

import java.sql.*;
import java.util.*;

import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class AthenaDataProvider extends JdbcDataProvider {
    public AthenaDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Athena";
        descriptor.description = "Query Athena database";
        descriptor.connectionTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.SERVER));
            add(new Property(Property.INT_TYPE, DbCredentials.PORT));
            add(new Property(Property.STRING_TYPE, DbCredentials.DB, DbCredentials.DB_DESCRIPTION));
            add(new Property(Property.STRING_TYPE, DbCredentials.LOGIN));
            add(new Property(Property.STRING_TYPE, DbCredentials.PASSWORD, new Prop("password")));
            add(new Property(Property.STRING_TYPE, "s3_staging_dir", "The name of the staging bucket"));
        }};
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("com.amazonaws.athena.jdbc.AthenaDriver");
        return DriverManager.getConnection(getConnectionString(conn), conn.getLogin(), conn.getPassword());
    }

    public String getConnectionString(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        String uri = "jdbc:awsathena://" + conn.getServer() + port;

        // Add supported parameters
        String parameters = "";
        String[] supportedParameters = {"s3_staging_dir"};
        for (int n = 0; n < supportedParameters.length; n++) {
            String key = supportedParameters[n];
            if (conn.parameters.keySet().contains(key)) {
                uri += (((parameters.equals(""))) ? "?" : "&") + key + "=" + conn.parameters.get(key);
            }
        }

        return uri + parameters;
    }
}

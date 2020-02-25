package grok_connect.providers;

import java.sql.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class RedshiftDataProvider extends JdbcDataProvider {
    public RedshiftDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Redshift";
        descriptor.category = "Database";
        descriptor.description = "Query Redshift database";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("com.amazon.redshift.jdbc42.Driver");
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (conn.parameters.containsKey(DbCredentials.SSL) && (boolean)conn.parameters.get(DbCredentials.SSL)) {
            properties.setProperty("ssl", "true");
            properties.setProperty("sslfactory", "com.amazon.redshift.ssl.NonValidatingFactory");
        }
        return DriverManager.getConnection(getConnectionString(conn), properties);
    }

    public String getConnectionString(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:redshift://" + conn.getServer() + port + "/" + conn.getDb();
    }
}
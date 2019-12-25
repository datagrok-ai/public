package grok_connect.providers;

import java.sql.*;
import grok_connect.connectors_info.*;


public class RedshiftDataProvider extends JdbcDataProvider {
    public RedshiftDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Redshift";
        descriptor.category = "Database";
        descriptor.description = "Query Redshift database";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("com.amazon.redshift.jdbc42.Driver");
        return DriverManager.getConnection(getConnectionString(conn), conn.getLogin(), conn.getPassword());
    }

    public String getConnectionString(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:redshift://" + conn.getServer() + port + "/" + conn.getDb();
    }
}
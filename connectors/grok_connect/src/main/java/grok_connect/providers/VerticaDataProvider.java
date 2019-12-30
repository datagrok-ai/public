package grok_connect.providers;

import java.sql.*;
import grok_connect.connectors_info.*;


public class VerticaDataProvider extends JdbcDataProvider {
    public VerticaDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Vertica";
        descriptor.description = "Query Vertica database";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("com.vertica.jdbc.Driver");
        return DriverManager.getConnection(getConnectionString(conn), conn.credentials.getLogin(), conn.credentials.getPassword());
    }

    public String getConnectionString(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:vertica://" + conn.getServer() + port + "/" + conn.getDb();
    }
}

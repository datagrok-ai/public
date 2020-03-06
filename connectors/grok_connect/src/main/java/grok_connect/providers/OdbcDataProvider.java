package grok_connect.providers;

import java.sql.*;
import grok_connect.connectors_info.*;


// TODO Note that this is just template. This class can be used as connector layer for other databases.
public class OdbcDataProvider extends JdbcDataProvider {
    public OdbcDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "ODBC";
        descriptor.description = "Query database via ODBC";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("sun.jdbc.odbc.JdbcOdbcDriver");
        return DriverManager.getConnection(getConnectionString(conn), conn.credentials.getLogin(), conn.credentials.getPassword());
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:odbc://" + conn.getServer() + port + "/" + conn.getDb();
    }
}

package grok_connect.providers;

import java.sql.*;
import grok_connect.connectors_info.*;


public class FirebirdDataProvider extends JdbcDataProvider {
    public FirebirdDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Firebird";
        descriptor.description = "Query Firebird database";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("org.firebirdsql.jdbc.FBDriver");
        return DriverManager.getConnection(getConnectionString(conn), conn.getLogin(), conn.getPassword());
    }

    public String getConnectionString(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : "/" + conn.getPort();
        return "jdbc:firebirdsql:" + conn.getServer() + port + ":" + conn.getDb();
    }
}

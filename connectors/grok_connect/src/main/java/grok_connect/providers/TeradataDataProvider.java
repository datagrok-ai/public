package grok_connect.providers;

import java.sql.*;
import grok_connect.connectors_info.*;


public class TeradataDataProvider extends JdbcDataProvider {
    public TeradataDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Teradata";
        descriptor.description = "Query Teradata database";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("com.teradata.jdbc.TeraDriver");
        return DriverManager.getConnection(getConnectionString(conn), conn.getLogin(), conn.getPassword());
    }

    public String getConnectionString(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:teradata://" + conn.getServer() + port + "/" + conn.getDb();
    }
}

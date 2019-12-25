package grok_connect.providers;

import java.sql.*;
import grok_connect.connectors_info.*;


public class HiveDataProvider extends JdbcDataProvider {
    public HiveDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Hive";
        descriptor.description = "Query Hive database";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("org.apache.hive.jdbc.HiveDriver");
        return DriverManager.getConnection(getConnectionString(conn), conn.getLogin(), conn.getPassword());
    }

    public String getConnectionString(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:hive://" + conn.getServer() + port + "/" + conn.getDb();
    }
}

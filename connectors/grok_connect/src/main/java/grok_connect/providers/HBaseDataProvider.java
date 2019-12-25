package grok_connect.providers;

import java.sql.*;
import grok_connect.connectors_info.*;


public class HBaseDataProvider extends JdbcDataProvider {
    public HBaseDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "HBase";
        descriptor.description = "Query HBase database";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("org.apache.phoenix.queryserver.client.Driver");
        return DriverManager.getConnection(getConnectionString(conn), conn.getLogin(), conn.getPassword());
    }

    public String getConnectionString(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:phoenix:thin:url=http://" + conn.getServer() + port + ";serialization=PROTOBUF";
    }
}

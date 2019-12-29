package grok_connect.providers;

import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;


public class VirtuosoDataProvider extends JdbcDataProvider {
    public VirtuosoDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Virtuoso";
        descriptor.description = "Query Virtuoso database";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("virtuoso.jdbc4.Driver");
        return DriverManager.getConnection(getConnectionString(conn), conn.getLogin(), conn.getPassword());
    }

    public String getConnectionString(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:virtuoso://" + conn.getServer() + port + "/TIMEOUT=100/UID=" + conn.getLogin() + "/PWD=" + conn.getPassword() + "/";
    }
}

package grok_connect.providers;

import java.sql.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class VirtuosoDataProvider extends JdbcDataProvider {
    public VirtuosoDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Virtuoso";
        descriptor.description = "Query Virtuoso database";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("virtuoso.jdbc4.Driver");
        return DriverManager.getConnection(getConnectionString(conn), conn.credentials.getLogin(), conn.credentials.getPassword());
    }

    public String getConnectionString(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        boolean ssl = conn.parameters.containsKey(DbCredentials.SSL) && (boolean)conn.parameters.get(DbCredentials.SSL);
        return "jdbc:virtuoso://" + conn.getServer() + port +
                "/TIMEOUT=100/UID=" + conn.credentials.getLogin() +
                "/PWD=" + conn.credentials.getPassword() +
                (ssl ? "/SSL" : "") +
                "/";
    }
}

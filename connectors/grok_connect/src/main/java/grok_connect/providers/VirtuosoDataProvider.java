package grok_connect.providers;

import java.sql.*;
import java.util.ArrayList;

import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class VirtuosoDataProvider extends JdbcDataProvider {
    public VirtuosoDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Virtuoso";
        descriptor.description = "Query Virtuoso database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("virtuoso.jdbc4.Driver");
        String connString = getConnectionString(conn);
        connString = connString.endsWith("/") ? connString : connString + "/";
        connString += "UID=" + conn.credentials.getLogin() + "/PWD=" + conn.credentials.getPassword();
        return DriverManager.getConnection(connString, conn.credentials.getLogin(), conn.credentials.getPassword());
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:virtuoso://" + conn.getServer() + port +
                "/TIMEOUT=100" + (conn.ssl() ? "/SSL" : "") + "/";
    }
}

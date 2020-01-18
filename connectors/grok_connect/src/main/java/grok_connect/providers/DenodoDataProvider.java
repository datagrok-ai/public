package grok_connect.providers;

import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.utils.Property;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.ArrayList;


public class DenodoDataProvider extends JdbcDataProvider {
    public DenodoDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Denodo";
        descriptor.description = "Query Impala database";

        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }


    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        //Class.forName("com.denodo.vdp.jdbc.Driver");
        Class.forName("com.denodo.vdb.jdbcdriver.VDBJDBCDriver");
        String connStr = getConnectionString(conn);
        return DriverManager.getConnection(connStr, conn.credentials.getLogin(), conn.credentials.getPassword());
    }

    public String getConnectionString(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:vdb://" + conn.getServer() + port +  "/" + conn.getDb();
    }
}

package grok_connect.providers;

import java.sql.*;
import org.apache.log4j.*;
import grok_connect.connectors_info.*;


public class CassandraDataProvider extends JdbcDataProvider {
    public CassandraDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Cassandra";
        descriptor.description = "Query Cassandra database";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("com.github.cassandra.jdbc.CassandraDriver");

        Logger logger = Logger.getLogger("com.github.cassandra.jdbc");
        logger.setLevel(Level.ERROR);

        return DriverManager.getConnection(getConnectionString(conn), conn.getLogin(), conn.getPassword());
    }

    public String getConnectionString(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:c*://" + conn.getServer() + port + "/" + conn.getDb();
    }
}

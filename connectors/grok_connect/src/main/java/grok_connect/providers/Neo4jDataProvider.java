package grok_connect.providers;

import java.sql.*;
import java.util.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class Neo4jDataProvider extends JdbcDataProvider {
    public Neo4jDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Neo4j";
        descriptor.description = "Query Neo4j database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.aggregations = null;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("org.neo4j.jdbc.Driver");
        Class.forName("org.neo4j.jdbc.bolt.BoltDriver");
        Class.forName("org.neo4j.jdbc.boltrouting.BoltRoutingNeo4jDriver");
        Class.forName("org.neo4j.jdbc.http.HttpDriver");
        java.util.Properties properties = defaultConnectionProperties(conn);
        String connString = getConnectionString(conn);
        if (!conn.hasCustomConnectionString() && !conn.ssl())
            connString += "?nossl";
        return DriverManager.getConnection(connString, properties);

    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:neo4j:bolt://" + conn.getServer() + port + "/" + conn.getDb();
    }
}

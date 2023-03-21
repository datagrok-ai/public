package grok_connect.providers;

import java.util.ArrayList;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.utils.Property;
import grok_connect.utils.ProviderManager;

public class Neo4jDataProvider extends JdbcDataProvider {
    public Neo4jDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "org.neo4j.jdbc.Driver";

        descriptor = new DataSource();
        descriptor.type = "Neo4j";
        descriptor.description = "Query Neo4j database";
        descriptor.commentStart = "//";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.aggregations = null;
    }

    @Override
    public String getConnectionString(DataConnection conn) {
        String connString = super.getConnectionString(conn);
        if (!conn.hasCustomConnectionString() && !conn.ssl())
            connString += "?nossl";
        return connString;
    }

    @Override
    public void prepareProvider() throws ClassNotFoundException {
        super.prepareProvider();
        Class.forName("org.neo4j.jdbc.bolt.BoltDriver");
        Class.forName("org.neo4j.jdbc.boltrouting.BoltRoutingNeo4jDriver");
        Class.forName("org.neo4j.jdbc.http.HttpDriver");
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:neo4j:bolt://" + conn.getServer() + port + "/" + conn.getDb();
    }
}

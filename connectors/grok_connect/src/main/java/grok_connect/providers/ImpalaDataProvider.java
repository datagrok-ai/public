package grok_connect.providers;

import java.sql.*;
import java.util.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class ImpalaDataProvider extends JdbcDataProvider {
    public ImpalaDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Impala";
        descriptor.description = "Query Impala database";
        descriptor.connectionTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.SERVER));
            add(new Property(Property.STRING_TYPE, DbCredentials.SCHEMA));
            add(new Property(Property.INT_TYPE, DbCredentials.PORT));
            add(new Property(Property.STRING_TYPE, DbCredentials.DB, DbCredentials.DB_DESCRIPTION));
            add(new Property(Property.STRING_TYPE, DbCredentials.CONNECTION_STRING,
                    DbCredentials.CONNECTION_STRING_DESCRIPTION, new Prop("textarea")));
        }};
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("com.cloudera.impala.jdbc41.Driver");
        return DriverManager.getConnection(getConnectionString(conn), conn.credentials.getLogin(), conn.credentials.getPassword());
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();

        String schema = (String)conn.parameters.get(DbCredentials.SCHEMA);
        schema = schema == null ? "/default" : "/" + schema;

        return "jdbc:impala://" + conn.getServer() + port + schema + ";AuthMech=3;SSL=1;AllowSelfSignedCerts=1";
    }
}

package grok_connect.providers;

import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.utils.Property;
import jdk.nashorn.internal.objects.Global;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.ArrayList;


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
        }};

        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("com.cloudera.impala.jdbc41.Driver");

        String connStr = getConnectionString(conn);
        System.out.println(connStr);

        return DriverManager.getConnection(connStr, conn.credentials.getLogin(), conn.credentials.getPassword());
    }

    public String getConnectionString(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();

        String schema = (String)conn.parameters.get(DbCredentials.SCHEMA);
        schema = schema == null ? "/default" : "/" + schema;

        return "jdbc:impala://" + conn.getServer() + port + schema + ";AuthMech=3;UID=" + conn.credentials.getLogin() + ";PWD=" + conn.credentials.getPasswordUrlEncoded();
    }
}

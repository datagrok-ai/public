package grok_connect.providers;

import java.sql.*;
import java.util.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class RedshiftDataProvider extends JdbcDataProvider {
    public RedshiftDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Redshift";
        descriptor.category = "Database";
        descriptor.description = "Query Redshift database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("com.amazon.redshift.jdbc42.Driver");
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString() && conn.ssl()) {
            properties.setProperty("ssl", "true");
            properties.setProperty("sslfactory", "com.amazon.redshift.ssl.NonValidatingFactory");
        }
        return DriverManager.getConnection(getConnectionString(conn), properties);
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:redshift://" + conn.getServer() + port + "/" + conn.getDb();
    }

    public String getSchemasSql(String db) {
        return "SELECT DISTINCT table_schema FROM information_schema.columns ORDER BY table_schema";
    }

    public String getSchemaSql(String db, String schema, String table) {
        List<String> filters = new ArrayList<String>() {{
            add("table_schema = '" + ((schema != null) ? schema : "public") + "'");
        }};

        if (db != null && db.length() != 0)
            filters.add("table_catalog = '" + db + "'");

        if (table != null)
            filters.add("table_name = '" + table + "'");

        String whereClause = "WHERE " + String.join(" AND \n", filters);

        return "SELECT table_schema, table_name, column_name, data_type " +
                "FROM information_schema.columns " + whereClause +
                " ORDER BY table_name";
    }
}
package grok_connect.providers;

import java.sql.*;
import java.util.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class TeradataDataProvider extends JdbcDataProvider {
    public TeradataDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Teradata";
        descriptor.description = "Query Teradata database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("com.teradata.jdbc.TeraDriver");
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString() && conn.ssl())
            properties.setProperty("ENABLESSL", "true");
        return DriverManager.getConnection(getConnectionString(conn), properties);

    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:teradata://" + conn.getServer() + port + "/" + conn.getDb();
    }

    public String getSchemasSql(String db) {
        return "SELECT DISTINCT databaseName as table_schema FROM DBC.TablesV ORDER BY databaseName";
    }

    public String getSchemaSql(String db, String schema, String table) {
        List<String> filters = new ArrayList<>();

        if (db == null || db.length() == 0)
            db = schema;

        if (db != null && db.length() != 0)
            filters.add("databaseName = '" + db + "'");

        if (table != null)
            filters.add("tableName = '" + table + "'");

        String whereClause = filters.size() != 0 ? "WHERE " + String.join(" and \n", filters) : "";

        return "SELECT databaseName as table_schema, tableName as table_name, " +
                "columnName as column_name, ColumnType as data_type " +
                "FROM DBC.ColumnsV " + whereClause + " ORDER BY tableName";
    }
}

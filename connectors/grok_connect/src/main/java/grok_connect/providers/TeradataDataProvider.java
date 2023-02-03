package grok_connect.providers;

import java.sql.*;
import java.util.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;
import serialization.Types;


public class TeradataDataProvider extends JdbcDataProvider {
    public TeradataDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "com.teradata.jdbc.TeraDriver";

        descriptor = new DataSource();
        descriptor.type = "Teradata";
        descriptor.description = "Query Teradata database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.canBrowseSchema = true;
        descriptor.nameBrackets = "\"";

        descriptor.typesMap = new HashMap<String, String>() {{
            put("bo", serialization.Types.BLOB);
        }};
    }

    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString() && conn.ssl())
            properties.setProperty("ENABLESSL", "true");
        return properties;

    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:teradata://" + conn.getServer() + port + "/" + conn.getDb();
    }

    public String getSchemasSql(String db) {
        return "SELECT DISTINCT databaseName as table_schema FROM DBC.TablesV";
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

    public String limitToSql(String query, Integer limit) {
        return query + "sample " + limit.toString();
    }

    public String addBrackets(String name) {
        String brackets = descriptor.nameBrackets;
        return name.startsWith(brackets.substring(0, 1)) ? name :
                brackets.substring(0, 1) + name + brackets.substring(brackets.length() - 1, brackets.length());
    }
}

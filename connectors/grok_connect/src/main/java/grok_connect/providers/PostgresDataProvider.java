package grok_connect.providers;

import java.sql.*;
import java.util.*;
import serialization.Types;
import grok_connect.utils.*;
import grok_connect.table_query.*;
import grok_connect.connectors_info.*;


public class PostgresDataProvider extends JdbcDataProvider {
    public PostgresDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "org.postgresql.Driver";

        descriptor = new DataSource();
        descriptor.type = "Postgres";
        descriptor.description = "Query Postgres database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.nameBrackets = "\"";

        descriptor.canBrowseSchema = true;
        descriptor.defaultSchema = "public";
        descriptor.typesMap = new HashMap<String, String>() {{
            put("smallint", Types.INT);
            put("int", Types.INT);
            put("bigint", Types.BIG_INT);
            put("real", Types.FLOAT);
            put("double precision", Types.FLOAT);
            put("numeric", Types.FLOAT);
            put("#character.*", Types.STRING);
            put("#varchar.*", Types.STRING);
            put("text", Types.STRING);
        }};
        descriptor.aggregations.add(new AggrFunctionInfo(Stats.STDEV, "stddev(#)", Types.dataFrameNumericTypes));
    }

    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString() && conn.ssl()) {
            properties.setProperty("ssl", "true");
            properties.setProperty("sslfactory", "org.postgresql.ssl.NonValidatingFactory");
        }
        return properties;
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:postgresql://" + conn.getServer() + port + "/" + conn.getDb();
    }

    public String getSchemasSql(String db) {
        return "SELECT DISTINCT table_schema FROM information_schema.columns";
    }

    public String getSchemaSql(String db, String schema, String table)
    {
        List<String> filters = new ArrayList<String>() {{
            add("table_schema = '" + ((schema != null) ? schema : descriptor.defaultSchema) + "'");
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

    private boolean isPostgresNumeric(String typeName) {
        // We need next condition because be default Postgres sets precision and scale to null for numeric type.
        // And ResultSetMetaData.getScale() returns 0 if scale is null.
        return typeName.equalsIgnoreCase("numeric");
    }

    @Override
    protected boolean isInteger(int type, String typeName, int precision, int scale) {
        if (isPostgresNumeric(typeName)) return false;
        return super.isInteger(type, typeName, precision, scale);
    }

    @Override
    protected boolean isFloat(int type, String typeName, int precision, int scale) {
        return super.isFloat(type, typeName, precision, scale) || isPostgresNumeric(typeName);
    }
}

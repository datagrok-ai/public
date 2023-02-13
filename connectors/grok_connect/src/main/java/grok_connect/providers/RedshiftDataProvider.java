package grok_connect.providers;

import java.sql.*;
import java.util.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;
import serialization.Types;
import grok_connect.table_query.*;

public class RedshiftDataProvider extends JdbcDataProvider {
    public RedshiftDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "com.amazon.redshift.jdbc.Driver";

        descriptor = new DataSource();
        descriptor.type = "Redshift";
        descriptor.category = "Database";
        descriptor.description = "Query Redshift database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;

        // copy-pasted from PostgresDataProvider; perhaps it would make sense to subclass it
        descriptor.nameBrackets = "\"";
        descriptor.canBrowseSchema = true;
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
            properties.setProperty("sslfactory", "com.amazon.redshift.ssl.NonValidatingFactory");
        }
        return properties;
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:redshift://" + conn.getServer() + port + "/" + conn.getDb();
    }

    public String getSchemasSql(String db) {
        return "SELECT DISTINCT table_schema FROM information_schema.columns";
    }

    public String getSchemaSql(String db, String schema, String table) {
        List<String> filters = new ArrayList<>();

        if (schema != null)
            filters.add("table_schema = '" + schema + "'");

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
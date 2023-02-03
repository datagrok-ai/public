package grok_connect.providers;

import java.sql.*;
import java.util.*;
import serialization.Types;
import grok_connect.utils.*;
import grok_connect.table_query.*;
import grok_connect.connectors_info.*;


public class MySqlDataProvider extends JdbcDataProvider {
    public MySqlDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "com.mysql.cj.jdbc.Driver";

        descriptor = new DataSource();
        descriptor.type = "MySQL";
        descriptor.description = "Query MySQL database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.canBrowseSchema = true;
        descriptor.nameBrackets = "`";

        descriptor.typesMap = new HashMap<String, String>() {{
            put("bool", Types.BOOL);
            put("boolean", Types.BOOL);
            put("#bit.*", Types.INT);
            put("#tinyint.*", Types.INT);
            put("#smallint.*", Types.INT);
            put("#mediumint.*", Types.INT);
            put("#int.*", Types.INT);
            put("#integer.*", Types.INT);
            put("#bigint.*", Types.BIG_INT);
            put("#dec.*", Types.FLOAT);
            put("#decimal.*", Types.FLOAT);
            put("#float.*", Types.FLOAT);
            put("#double.*", Types.FLOAT);
            put("#double precision.*", Types.FLOAT);
            put("#char.*", Types.STRING);
            put("#varchar.*", Types.STRING);
            put("#text.*", Types.STRING);
            put("#tinytext.*", Types.STRING);
            put("date", Types.DATE_TIME);
            put("#datetime.*", Types.DATE_TIME);
            put("#timestamp.*", Types.DATE_TIME);
            put("#time.*", Types.DATE_TIME);
        }};
        descriptor.aggregations.add(new AggrFunctionInfo(Stats.STDEV, "std(#)", Types.dataFrameNumericTypes));
    }

    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString()) {
            properties.setProperty("zeroDateTimeBehavior", "convertToNull");
            if (conn.ssl()) {
                properties.setProperty("useSSL", "true");
                properties.setProperty("verifyServerCertificate", "false");
            }
        }
        return properties;
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:mysql://" + conn.getServer() + port + "/" + conn.getDb();
    }

    public String getSchemasSql(String db) {
        return "SELECT DISTINCT table_schema FROM information_schema.columns";
    }

    public String getSchemaSql(String db, String schema, String table) {
        List<String> filters = new ArrayList<>();

        if (db == null || db.length() == 0)
            db = schema;

        if (db != null && db.length() != 0)
            filters.add("table_schema = '" + db + "'");

        if (table != null)
            filters.add("(table_name = '" + table + "')");

        String whereClause = filters.size() != 0 ? "WHERE " + String.join(" AND \n", filters) : "";

        return "SELECT table_schema, table_name, column_name, data_type " +
                "FROM information_schema.columns " + whereClause + " ORDER BY table_name";
    }
}

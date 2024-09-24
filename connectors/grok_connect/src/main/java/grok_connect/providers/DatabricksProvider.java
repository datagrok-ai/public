package grok_connect.providers;

import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.Prop;
import grok_connect.utils.Property;
import grok_connect.utils.QueryCancelledByUser;
import serialization.DataFrame;
import serialization.IntColumn;
import serialization.StringColumn;
import serialization.Types;
import java.sql.Connection;
import java.sql.DatabaseMetaData;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.*;

public class DatabricksProvider extends JdbcDataProvider {
    public DatabricksProvider() {
        driverClassName = "com.databricks.client.jdbc.Driver";

        descriptor = new DataSource();
        descriptor.type = "Databricks";
        descriptor.description = "Query Databricks warehouse";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.UID,
                    "Set the UID property to an appropriate user name for accessing the Spark server"));
            add(new Property(Property.STRING_TYPE, DbCredentials.PWD, "Set the PWD property to the password corresponding to the user name you\n" +
                    "provided. If UID is either set to token or left at the default value (which is also\n" +
                    "token), then set the PWD property to the Databricks token value that you want to use for authentication.",
                    new Prop("password")));
        }};
        descriptor.nameBrackets = "`";
        descriptor.defaultSchema = "default";
        descriptor.canBrowseSchema = true;
        descriptor.typesMap = new HashMap<String, String>() {{
            put("tinyint", Types.INT);
            put("smallint", Types.INT);
            put("integer", Types.INT);
            put("int", Types.INT);
            put("bigint", Types.BIG_INT);
            put("real", Types.FLOAT);
            put("double", Types.FLOAT);
            put("#decimal.*", Types.FLOAT);
            put("boolean", Types.BOOL);
            put("#char.*", Types.STRING);
            put("#varchar.*", Types.STRING);
            put("string", Types.STRING);
            put("text", Types.STRING);
            put("date", Types.DATE_TIME);
            put("#timestamp.*", Types.DATE_TIME);
            put("#varbinary.*", Types.BLOB);
        }};
    }

    @Override
    public Properties getProperties(DataConnection conn) {
        Properties properties = new Properties();
        setIfNotNull(properties, "user", (String) conn.credentials.parameters.get(DbCredentials.UID));
        setIfNotNull(properties, "password", (String) conn.credentials.parameters.get(DbCredentials.PWD));
        if (!conn.hasCustomConnectionString() && conn.ssl())
            properties.setProperty("SSL", "1");
        return properties;
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:databricks://" + conn.getServer() + port +  "/" + conn.getDb();
    }

    @Override
    public DataFrame getSchemas(DataConnection connection) throws QueryCancelledByUser, GrokConnectException {
        try (Connection dbConnection = getConnection(connection);
             ResultSet catalogs = dbConnection.getMetaData().getCatalogs()) {
            DataFrame result = DataFrame.fromColumns(new StringColumn("table_schema"));
            while (catalogs .next())
                result.addRow(catalogs.getString(1));
            return result;
        } catch (SQLException e) {
            throw new GrokConnectException(e);
        }
    }

    @Override
    public DataFrame getSchema(DataConnection connection, String schema, String table) throws
            QueryCancelledByUser, GrokConnectException {
        try (Connection dbConnection = getConnection(connection)) {
            IntColumn isView = new IntColumn("is_view");
            DataFrame result = DataFrame.fromColumns(new StringColumn("table_schema"),
                    new StringColumn("table_name"), new StringColumn("column_name"),
                    new StringColumn("data_type"), isView);
            Map<String, Integer> tableIndexes = new HashMap<>();
            DatabaseMetaData metaData = dbConnection.getMetaData();
            try (ResultSet columns = metaData.getColumns(schema, null, table, null)) {
                while (columns.next()) {
                    tableIndexes.put(columns.getString(3), result.rowCount);
                    result.addRow(columns.getString(1), columns.getString(3),
                            columns.getString(4), columns.getString(6), null);
                }
            }

            try (ResultSet tables = metaData.getTables(schema, null, table, null)) {
                while (tables.next()) {
                    Integer index = tableIndexes.get(tables.getString("TABLE_NAME"));
                    isView.set(index,  "VIEW".equals(tables.getString("TABLE_TYPE")) ? 1 : 0);
                }
            }
            return result;
        } catch (SQLException e) {
            throw new GrokConnectException(e);
        }
    }

    @Override
    public DataFrame getForeignKeys(DataConnection conn, String schema) throws GrokConnectException {
        try (Connection connection = getConnection(conn)) {
            DatabaseMetaData meta = connection.getMetaData();
            List<String> tables = new ArrayList<>();
            try (ResultSet tablesRs = meta.getTables(conn.getDb(), null, null, new String[]{"TABLE", "VIEW"})) {
                while (tablesRs.next())
                    tables.add(tablesRs.getString("TABLE_NAME"));
            }

            DataFrame result = DataFrame.fromColumns(new StringColumn("table_schema"),
                    new StringColumn("constraint_name"), new StringColumn("table_name"),
                    new StringColumn("column_name"), new StringColumn("foreign_table_name"), new StringColumn("foreign_column_name"));
            if (!tables.isEmpty()) {
                for (String t : tables)
                    try (ResultSet info = meta.getExportedKeys(conn.getDb(), null, t)) {
                        while(info.next())
                            result.addRow(info.getString("FKTABLE_SCHEM"), info.getString("FK_NAME"),
                                    info.getString("FKTABLE_NAME"), info.getString("FKCOLUMN_NAME"),
                                    info.getString("PKTABLE_NAME"), info.getString("PKCOLUMN_NAME"));
                    }
            }
            return result;
        } catch (SQLException e) {
            throw new GrokConnectException(e);
        }
    }
}

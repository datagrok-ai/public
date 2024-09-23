package grok_connect.providers;

import java.sql.Connection;
import java.sql.DatabaseMetaData;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;
import serialization.DataFrame;
import serialization.IntColumn;
import serialization.StringColumn;
import serialization.Types;

public class DenodoDataProvider extends JdbcDataProvider {
    public DenodoDataProvider() {
        driverClassName = "com.denodo.vdp.jdbc.Driver";

        descriptor = new DataSource();
        descriptor.type = "Denodo";
        descriptor.description = "Query Denodo database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.nameBrackets = "\"";
        descriptor.canBrowseSchema = true;
        descriptor.typesMap = new HashMap<String, String>() {{
            put("integer", Types.INT);
            put("smallint", Types.INT);
            put("tinyiny", Types.INT);
            put("bigint", Types.BIG_INT);
            put("boolean", Types.BOOL);
            put("bit", Types.BOOL);
            put("varchar", Types.STRING);
            put("text", Types.STRING);
            put("decimal", Types.FLOAT);
            put("float", Types.FLOAT);
            put("double", Types.FLOAT);
            put("xml", Types.OBJECT);
            put("time", Types.DATE_TIME);
            put("timestamp", Types.DATE_TIME);
            put("date", Types.DATE_TIME);
            put("blob", serialization.Types.BLOB);
            put("clob", serialization.Types.BLOB);
        }};
    }

    @Override
    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString() && conn.ssl())
            properties.setProperty("ssl", "true");
        return properties;
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:denodo://" + conn.getServer() + port +  "/" + conn.getDb();
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
}

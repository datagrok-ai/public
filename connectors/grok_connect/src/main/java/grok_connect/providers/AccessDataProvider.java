package grok_connect.providers;

import java.io.File;
import java.sql.Connection;
import java.sql.DatabaseMetaData;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.table_query.AggrFunctionInfo;
import grok_connect.table_query.GroupAggregation;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.Prop;
import grok_connect.utils.Property;
import serialization.DataFrame;
import serialization.StringColumn;
import serialization.Types;

public class AccessDataProvider extends JdbcDataProvider {
    public AccessDataProvider() {
        driverClassName = "net.ucanaccess.jdbc.UcanaccessDriver";
        descriptor = new DataSource();
        descriptor.type = "Access";
        descriptor.description = "Query Access database";
        descriptor.canBrowseSchema = true;
        descriptor.typesMap = new HashMap<String, String>() {{
            put("integer", serialization.Types.INT);
            put("smallint", serialization.Types.INT);
            put("boolean", Types.BOOL);
            put("decimal", serialization.Types.FLOAT);
            put("numeric", serialization.Types.FLOAT);
            put("double", serialization.Types.FLOAT);
            put("varchar", serialization.Types.STRING);
            put("timestamp", Types.DATE_TIME);
            put("blob", serialization.Types.BLOB);
        }};
        descriptor.connectionTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.CONNECTION_STRING,
                    DbCredentials.CONNECTION_STRING_DESCRIPTION, new Prop("textarea")));
            add(new Property(Property.STRING_TYPE, DbCredentials.DB, DbCredentials.DB_DESCRIPTION));
        }};
        descriptor.credentialsTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.LOGIN));
            add(new Property(Property.STRING_TYPE, DbCredentials.PASSWORD, new Prop("password")));
        }};
    }

    @Override
    public DataFrame getSchemas(DataConnection connection) {
        StringColumn column = new StringColumn(new String[]{""});
        column.name = "TABLE_SCHEMA";
        DataFrame dataFrame = new DataFrame();
        dataFrame.addColumn(column);
        return dataFrame;
    }

    @Override
    public DataFrame getSchema(DataConnection connection, String schema, String table) throws GrokConnectException {
        try (Connection dbConnection = getConnection(connection);
             ResultSet columns = dbConnection.getMetaData().getColumns(null, schema, table, null)) {
            DataFrame result = DataFrame.fromColumns(new StringColumn("table_schema"),
                    new StringColumn("table_name"), new StringColumn("column_name"),
                    new StringColumn("data_type"));
            while (columns.next())
                result.addRow(columns.getString(2), columns.getString(3),
                        columns.getString(4), columns.getString(6));
            return result;
        } catch (SQLException e) {
            throw new GrokConnectException(e);
        }
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        return String.format("jdbc:ucanaccess://%s;memory=false", new File(conn.getDb()).getAbsolutePath());
    }

    @Override
    public DataFrame getForeignKeys(DataConnection conn, String schema) throws GrokConnectException {
        try (Connection connection = getConnection(conn)) {
            DatabaseMetaData meta = connection.getMetaData();
            List<String> tables = new ArrayList<>();
            try (ResultSet tablesRs = meta.getTables(conn.getDb(), null, null, null)) {
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

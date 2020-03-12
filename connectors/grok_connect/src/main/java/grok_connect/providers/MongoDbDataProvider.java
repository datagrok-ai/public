package grok_connect.providers;

import serialization.Types;
import org.bson.*;
import java.sql.*;
import java.util.*;
import serialization.*;
import grok_connect.connectors_info.*;


public class MongoDbDataProvider extends JdbcDataProvider {
    public MongoDbDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "MongoDB";
        descriptor.description = "Query MongoDB database";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("com.dbschema.MongoJdbcDriver");
        return DriverManager.getConnection(getConnectionString(conn), conn.credentials.getLogin(), conn.credentials.getPassword());
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:mongodb://" + conn.getServer() + port + "/" + conn.getDb();
    }

    @SuppressWarnings("unchecked")
    public DataFrame execute(FuncCall queryRun)
            throws ClassNotFoundException, SQLException {

        DataQuery dataQuery = queryRun.func;
        Connection connection = getConnection(dataQuery.connection);
        Statement statement = connection.createStatement();
        ResultSet resultSet = statement.executeQuery(dataQuery.query);

        int rowCount = 0;
        int columnCount = 0;
        List<Column> columns = new ArrayList<>(columnCount);

        while (resultSet.next()) {

            Document collection = (Document)resultSet.getObject(1);
            Set<String> columnNames = collection.keySet();

            if (rowCount++ <= 0) {
                // Detect all columns form first collection
                for (String columnName : columnNames) {
                    Column column;
                    Object value = collection.get(columnName);
                    if ((value instanceof Byte) || (value instanceof Short) || (value instanceof Integer))
                        column = new IntColumn();
                    else if (value instanceof Float)
                        column = new FloatColumn();
                    else if (value instanceof Double) {
                        column = new FloatColumn();
                        value = new Float((Double)value);
                    } else if ((value instanceof Boolean))
                        column = new BoolColumn();
                    else {
                        column = new StringColumn();
                        value = value.toString();
                    }

                    column.add(value);
                    column.name = columnName;
                    columns.add(column);
                }
            } else {
                // Process all next collections
                for (Column column : columns) {
                    Object value = collection.get(column.name);
                    String colType = column.getType();
                    switch (colType) {
                        case Types.FLOAT: {
                            if (value instanceof Double)
                                column.add(new Float((Double) value));
                            else
                                column.add(value);
                            break;
                        }
                        case Types.INT:
                        case Types.BOOL:
                            column.add(value);
                            break;
                        default:
                            column.add((value != null) ? value.toString() : "");
                    }
                }
            }
        }

        connection.close();

        DataFrame dataFrame = new DataFrame();
        dataFrame.addColumns(columns);

        return dataFrame;
    }
}

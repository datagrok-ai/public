package grok_connect.providers;

import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataQuery;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncCall;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.ProviderManager;
import grok_connect.utils.SchemeInfo;
import org.bson.Document;
import serialization.BoolColumn;
import serialization.Column;
import serialization.DataFrame;
import serialization.FloatColumn;
import serialization.IntColumn;
import serialization.StringColumn;
import serialization.Types;
import java.sql.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

public class MongoDbDataProvider extends JdbcDataProvider {
    public MongoDbDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "com.dbschema.MongoJdbcDriver";

        descriptor = new DataSource();
        descriptor.type = "MongoDB";
        descriptor.description = "Query MongoDB database";
        descriptor.connectionTemplate = DbCredentials.dbConnectionTemplate;
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.aggregations = null;
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:mongodb://" + conn.getServer() + port + "/" + conn.getDb();
    }

    @Override
    @SuppressWarnings("unchecked")
    public DataFrame execute(FuncCall queryRun)
            throws ClassNotFoundException, SQLException, GrokConnectException {

        DataQuery dataQuery = queryRun.func;
        Connection connection = getConnection(dataQuery.connection);
        ResultSet resultSet = getResultSet(queryRun, connection);
        List<Column> columns = getTypedColumns(resultSet);
        DataFrame dataFrame = getResultSetSubDf(queryRun, resultSet, columns,
                new ArrayList<>(), new ArrayList<>(), 0);
        return dataFrame;
    }

    @Override
    public ResultSet getResultSet(FuncCall queryRun, Connection connection) {
        try {
            PreparedStatement statement = connection.prepareStatement(queryRun.func.query);
            return statement.executeQuery();
        } catch (SQLException e) {
            throw new RuntimeException("Something went wrong when getting resultSet", e);
        }
    }

    @Override
    public SchemeInfo resultSetScheme(FuncCall queryRun, ResultSet resultSet) {
        return new SchemeInfo(getTypedColumns(resultSet), new ArrayList<>(), new ArrayList<>());
    }

    @Override
    @SuppressWarnings("unchecked")
    public DataFrame getResultSetSubDf(FuncCall queryRun, ResultSet resultSet, List<Column> columns,
                                    List<Boolean> supportedType,List<Boolean> initColumn, int maxIterations) throws SQLException {
        do {
            Object object = resultSet.getObject(1);
            if (object instanceof String) {
                columns.get(0).add(object.toString());
                continue;
            }
            Document collection = (Document) object;
            if (collection == null) {
                continue;
            }
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
        } while (resultSet.next());
        resultSet.close();
        DataFrame dataFrame = new DataFrame();
        dataFrame.addColumns(columns);
        return dataFrame;
    }

    private List<Column> getTypedColumns(ResultSet resultSet) {
        List<Column> columns = new ArrayList<>();
        Object object;
        String label;
        if (resultSet == null) {
            return columns;
        }
        try {
            resultSet.next();
            object = resultSet.getObject(1);
            label = resultSet.getMetaData().getColumnLabel(1);
        } catch (SQLException e) {
            throw new RuntimeException("Something went wrong when retrieving meta data of resultSet", e);
        }
        if (object == null) {
            return columns;
        } else if (object instanceof String) {
            StringColumn column = new StringColumn();
            column.name = label;
            columns.add(column);
        } else {
            Document collection = (Document) object;
            Set<String> columnNames = collection.keySet();
            for (String columnName : columnNames) {
                Column column;
                Object value = collection.get(columnName);
                if ((value instanceof Byte) || (value instanceof Short) || (value instanceof Integer)) {
                    column = new IntColumn();
                } else if (value instanceof Float)
                    column = new FloatColumn();
                else if (value instanceof Double) {
                    column = new FloatColumn();
                } else if ((value instanceof Boolean))
                    column = new BoolColumn();
                else {
                    column = new StringColumn();
                }
                column.name = columnName;
                columns.add(column);
            }
        }
        return columns;
    }
}

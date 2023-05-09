package grok_connect.providers;

import grok_connect.column.ColumnManager;
import grok_connect.column.complex.Neo4jComplexColumnManager;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataQuery;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncCall;
import grok_connect.resultset.DefaultResultSetManager;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.SchemeInfo;
import serialization.Column;
import serialization.DataFrame;
import serialization.Types;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class MongoDbDataProvider extends JdbcDataProvider {
    private static final int OBJECT_INDEX = 1;

    public MongoDbDataProvider() {
        initResultSetManager();
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
        while (resultSet.next()) {
            Object object = resultSet.getObject(OBJECT_INDEX);
            if (object instanceof String) {
                columns.get(0).add(object.toString());
                continue;
            }
            Column column = columns.get(0);
            column.add(resultSetManager.convert(object, 2000,
                    resultSet.getMetaData().getColumnTypeName(OBJECT_INDEX), 0, 0, column.name));
        };
        resultSet.close();
        DataFrame dataFrame = new DataFrame();
        dataFrame.addColumns(columns);
        return dataFrame;
    }

    private List<Column> getTypedColumns(ResultSet resultSet) {
        List<Column> columns = new ArrayList<>();
        String label;
        String typeName;
        int type;
        if (resultSet == null) {
            return columns;
        }
        try {
            ResultSetMetaData metaData = resultSet.getMetaData();
            label = metaData.getColumnLabel(OBJECT_INDEX);
            type = metaData.getColumnType(OBJECT_INDEX);
            typeName = metaData.getColumnTypeName(OBJECT_INDEX);
            Column column = resultSetManager.getColumn(type, typeName, 0, 0);
            column.name = label;
            columns.add(column);
            return columns;
        } catch (SQLException e) {
            throw new RuntimeException("Something went wrong when retrieving meta data of resultSet", e);
        }
    }

    private void initResultSetManager() {
        Map<String, ColumnManager<?>> defaultManagersMap = DefaultResultSetManager.getDefaultManagersMap();
        defaultManagersMap.put(Types.COLUMN_LIST, new Neo4jComplexColumnManager());
        resultSetManager = DefaultResultSetManager.fromManagersMap(defaultManagersMap);
    }
}

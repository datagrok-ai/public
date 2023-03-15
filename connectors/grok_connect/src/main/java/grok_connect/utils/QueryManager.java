package grok_connect.utils;

import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.*;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.JdbcDataProvider;
import serialization.*;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import org.apache.log4j.Logger;
import org.joda.time.DateTime;

public class QueryManager {
    JdbcDataProvider provider;

    SchemeInfo schemeInfo;

    public ResultSet resultSet;
    public FuncCall query;
    Connection connection;

    static Gson gson = new GsonBuilder()
        .registerTypeAdapter(Property.class, new PropertyAdapter())
        .create();

    private static ProviderManager providerManager = new ProviderManager(Logger.getLogger(GrokConnect.class.getName()));

    public QueryManager(FuncCall query, JdbcDataProvider provider) {
        this.query = query;
        this.provider = provider;
    }

    public QueryManager(String message) {
        query = gson.fromJson(message, FuncCall.class);
        query.log = "";
        query.setParamValues();
        query.afterDeserialization();
        System.out.println(query.func.query);

        // DateTime startTime = DateTime.now();
        provider = providerManager.getByName(query.func.connection.dataSource);
    }

    public void getResultSet() throws ClassNotFoundException, GrokConnectException, QueryCancelledByUser, SQLException {
        connection = provider.getConnection(query.func.connection);
        resultSet = provider.getResultSet(query, connection);
        System.out.println("resultset got");
    }

    public void initScheme() throws QueryCancelledByUser, SQLException {
        System.out.println("scheme not init");
        schemeInfo = provider.resultSetScheme(query, resultSet);
        System.out.println("scheme init");
    }

    public DataFrame getSubDF(int maxIterations) throws IOException, SQLException, QueryCancelledByUser {
        List<Column> columns = schemeInfo.columns;
        for (int i = 0; i < columns.size(); i++) {
            columns.get(i).empty();
        }

        DataFrame df = provider.getResultSetSubDf(query, resultSet, columns, schemeInfo.supportedType, schemeInfo.initColumn, maxIterations);
        return df;
    }

    public void closeConnection() throws SQLException {
        if (connection != null && !connection.isClosed()) {
            connection.commit();
            connection.close();
        }
    }
}

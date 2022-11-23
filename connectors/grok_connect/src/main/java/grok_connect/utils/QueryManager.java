package grok_connect.utils;

import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.*;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.JdbcDataProvider;
import serialization.*;;

public class QueryManager {
    JdbcDataProvider provider;

    SchemeInfo schemeInfo;

    ResultSet resultSet;
    FuncCall query;


    public QueryManager(FuncCall query, JdbcDataProvider provider) {
        this.query = query;
        this.provider = provider;
    }

    public void getResultSet() throws ClassNotFoundException, GrokConnectException, QueryCancelledByUser, SQLException {
        resultSet = provider.getResultSet(query);
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

        return provider.getResultSetSubDf(query, resultSet, columns, schemeInfo.supportedType, schemeInfo.initColumn, maxIterations);
    }
    

}

package grok_connect.providers;

import java.sql.*;
import java.util.*;

import grok_connect.resultset.ResultSetManager;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;
import serialization.Types;


public class ImpalaDataProvider extends JdbcDataProvider {
    public ImpalaDataProvider(ResultSetManager resultSetManager, ProviderManager providerManager) {
        super(resultSetManager, providerManager);
        driverClassName = "com.cloudera.impala.jdbc41.Driver";

        descriptor = new DataSource();
        descriptor.type = "Impala";
        descriptor.description = "Query Impala database";
        descriptor.connectionTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.SERVER));
            add(new Property(Property.STRING_TYPE, DbCredentials.SCHEMA));
            add(new Property(Property.INT_TYPE, DbCredentials.PORT, new Prop()));
            add(new Property(Property.STRING_TYPE, DbCredentials.DB, DbCredentials.DB_DESCRIPTION));
            add(new Property(Property.STRING_TYPE, DbCredentials.CONNECTION_STRING,
                    DbCredentials.CONNECTION_STRING_DESCRIPTION, new Prop("textarea")));
            add(new Property(Property.BOOL_TYPE, DbCredentials.CACHE_SCHEMA));
            add(new Property(Property.BOOL_TYPE, DbCredentials.CACHE_RESULTS));
            add(new Property(Property.STRING_TYPE, DbCredentials.CACHE_INVALIDATE_SCHEDULE));
        }};
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    protected void appendQueryParam(DataQuery dataQuery, String paramName, StringBuilder queryBuffer) {
        //if list -- append all items
        FuncParam param = dataQuery.getParam(paramName);
        if (param.propertyType.equals(Types.LIST)) {
            if (param.value == null) {
                queryBuffer.append("?");
                return;
            }
            @SuppressWarnings (value="unchecked")
            ArrayList<Object> lst = (ArrayList<Object>)param.value;
            int size = lst.size();
            if (size == 0) {
                queryBuffer.append("?");
                return;
            }
            for (int i = 0; i < size; i++) {
                queryBuffer.append("?");
                if (i < size - 1)
                    queryBuffer.append(",");
            }
        } else {
            queryBuffer.append("?");
        }
    }

    protected int setArrayParamValue(PreparedStatement statement, int n, FuncParam param) throws SQLException {
        //iterate ist and add all the parameters
        @SuppressWarnings (value="unchecked")
        ArrayList<Object> lst = (ArrayList<Object>)param.value;
        if (lst == null || lst.size() == 0) {
            statement.setObject(n, null);
            return 0;
        }
        for (int i = 0; i < lst.size(); i++) {
            System.out.println(n + i);
            statement.setObject(n + i, lst.get(i));
        }
        return lst.size() - 1;
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();

        String schema = (String)conn.parameters.get(DbCredentials.SCHEMA);
        schema = schema == null ? "/default" : "/" + schema;

        return "jdbc:impala://" + conn.getServer() + port + schema + ";AuthMech=3;SSL=1;AllowSelfSignedCerts=1";
    }
}

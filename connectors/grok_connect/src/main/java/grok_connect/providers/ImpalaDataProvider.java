package grok_connect.providers;

import java.sql.*;
import java.util.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;
import serialization.Types;


public class ImpalaDataProvider extends JdbcDataProvider {
    public ImpalaDataProvider(ProviderManager providerManager) {
        super(providerManager);
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
        }};
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName(driverClassName);
        return CustomDriverManager.getConnection(getConnectionString(conn), conn.credentials.getLogin(), conn.credentials.getPassword(), driverClassName);
    }

    protected void appendQueryParam(DataQuery dataQuery, String paramName, StringBuilder queryBuffer) {
        //if list -- append all items
        FuncParam param = dataQuery.getParam(paramName);
        if (param.propertyType.equals(Types.LIST)) {
            @SuppressWarnings (value="unchecked")
            ArrayList<Object> lst = (ArrayList<Object>)param.value;
            int size = lst.size();
            for (int i = 0; i < size; i++) {
                queryBuffer.append("?");
                if (i < size - 1)
                    queryBuffer.append(",");
            }
        } else {
            queryBuffer.append("?");
        }
    }

    protected void setArrayParamValue(PreparedStatement statement, int n, FuncParam param) throws SQLException {
        //iterate ist and add all the parameters
        @SuppressWarnings (value="unchecked")
        ArrayList<Object> lst = (ArrayList<Object>)param.value;
        for (int i = 0; i < lst.size(); i++) {
            statement.setObject(n + 1 + i, lst.get(i));
        }
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();

        String schema = (String)conn.parameters.get(DbCredentials.SCHEMA);
        schema = schema == null ? "/default" : "/" + schema;

        return "jdbc:impala://" + conn.getServer() + port + schema + ";AuthMech=3;SSL=1;AllowSelfSignedCerts=1";
    }
}

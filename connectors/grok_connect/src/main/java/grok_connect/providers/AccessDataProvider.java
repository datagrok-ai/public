package grok_connect.providers;

import java.sql.*;
import java.util.*;

import serialization.Types;
import grok_connect.utils.*;
import grok_connect.table_query.*;
import grok_connect.connectors_info.*;


public class AccessDataProvider extends JdbcDataProvider {
    public AccessDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Access";
        descriptor.description = "Query Access database";
        descriptor.connectionTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.CONNECTION_STRING,
                    DbCredentials.CONNECTION_STRING_DESCRIPTION, new Prop("textarea")));
            add(new Property(Property.STRING_TYPE, DbCredentials.DB, DbCredentials.DB_DESCRIPTION));
        }};
        descriptor.credentialsTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.LOGIN));
            add(new Property(Property.STRING_TYPE, DbCredentials.PASSWORD, new Prop("password")));
        }};
        descriptor.aggregations.add(new AggrFunctionInfo(Stats.STDEV, "stdev(#)", Types.dataFrameNumericTypes));
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("net.ucanaccess.jdbc.UcanaccessDriver");
        return DriverManager.getConnection(getConnectionString(conn), conn.credentials.getLogin(),
                conn.credentials.getPassword());
    }

    public String getConnectionStringImpl(DataConnection conn) {
        return "jdbc:ucanaccess://" + conn.getDb();
    }
}

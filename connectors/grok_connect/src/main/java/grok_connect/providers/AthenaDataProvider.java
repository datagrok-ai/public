package grok_connect.providers;

import java.sql.*;
import java.util.*;

import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class AthenaDataProvider extends JdbcDataProvider {
    public AthenaDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Athena";
        descriptor.description = "Query Athena database";
        descriptor.connectionTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.SERVER));
            add(new Property(Property.INT_TYPE, DbCredentials.PORT, new Prop()));
            add(new Property(Property.STRING_TYPE, DbCredentials.DB, DbCredentials.DB_DESCRIPTION));
            add(new Property(Property.STRING_TYPE, DbCredentials.CONNECTION_STRING,
                    DbCredentials.CONNECTION_STRING_DESCRIPTION, new Prop("textarea")));
            add(new Property(Property.STRING_TYPE, "S3OutputLocation",
                    "The path of the Amazon S3 location where you want to store query results"));
        }};
        descriptor.credentialsTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, "AccessKey"));
            add(new Property(Property.STRING_TYPE, "SecretKey", new Prop("password")));
        }};
    }

    public boolean isParametrized() {
        return false;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("com.simba.athena.jdbc.Driver");
        String connString = getConnectionString(conn);
        connString = connString.endsWith(";") ? connString : connString + ";";
        connString += "User=" + conn.credentials.parameters.get("AccessKey") + ";" +
                "Password=" + conn.credentials.parameters.get("SecretKey");
        return DriverManager.getConnection(connString);
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:awsathena://" + conn.getServer() + port + ";" +
                "S3OutputLocation=" + conn.parameters.get("S3OutputLocation");
    }

    public String castParamValueToSqlDateTime(FuncParam param) {
        return "from_iso8601_timestamp('" + param.value.toString() + "')";
    }
}

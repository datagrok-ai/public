package grok_connect.providers;

import java.sql.*;
import java.util.*;

import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class AthenaDataProvider extends JdbcDataProvider {
    public AthenaDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "com.simba.athena.jdbc41.Driver";

        Property encode = new Property(Property.STRING_TYPE, "s3OutputEncOption",
                "The encryption protocol that the driver uses to encrypt your query results before storing them on Amazon S3");
        encode.choices = new ArrayList<String>() {{ add("SSE_S3"); add("SSE_KMS"); add("CSE_KMS"); }};

        descriptor = new DataSource();
        descriptor.type = "Athena";
        descriptor.description = "Query Athena database";
        descriptor.connectionTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.SERVER));
            add(new Property(Property.INT_TYPE, DbCredentials.PORT, new Prop()));
            add(new Property(Property.STRING_TYPE, DbCredentials.DB, DbCredentials.DB_DESCRIPTION));
            add(new Property(Property.STRING_TYPE, "s3OutputLocation",
                    "The path of the Amazon S3 location where you want to store query results"));
            add(encode);
            add(new Property(Property.STRING_TYPE, DbCredentials.CONNECTION_STRING,
                    DbCredentials.CONNECTION_STRING_DESCRIPTION, new Prop("textarea")));
        }};
        descriptor.credentialsTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, "accessKey"));
            add(new Property(Property.STRING_TYPE, "secretKey", new Prop("password")));
        }};
    }

    public boolean autoInterpolation() {
        return false;
    }


    public String getConnectionString(DataConnection conn) {
        String connString = super.getConnectionString(conn);
        connString = connString.endsWith(";") ? connString : connString + ";";
        String accessKey = "";
        String secretKey = "";
        if (conn.credentials != null) {
            if (conn.credentials.parameters.get("accessKey") != null)
                accessKey = conn.credentials.parameters.get("accessKey").toString();
            if (conn.credentials.parameters.get("secretKey") != null)
                secretKey = conn.credentials.parameters.get("secretKey").toString();
        }
        connString += "User=" + accessKey + ";" +
                "Password=" + secretKey;
        return connString;
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        String encode = (String)conn.parameters.get("s3OutputEncOption");
        return "jdbc:awsathena://" + conn.getServer() + port +
                ";S3OutputLocation=" + conn.parameters.get("s3OutputLocation") +
                (encode == null ? "" : ";S3OutputEncOption=" + encode);
    }

    public String castParamValueToSqlDateTime(FuncParam param) {
        return "from_iso8601_timestamp('" + param.value.toString() + "')";
    }
}

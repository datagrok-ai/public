package grok_connect.providers;

import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncParam;
import grok_connect.utils.*;

import java.sql.Connection;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Properties;

public class NeptuneDataProvider extends JdbcDataProvider {
    public NeptuneDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "software.aws.neptune.NeptuneDriver";
        descriptor = new DataSource();
        descriptor.type = "Neptune";
        descriptor.description = "Amazon Neptune is a fully managed graph database service";
        descriptor.connectionTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.SERVER));
            add(new Property(Property.INT_TYPE, DbCredentials.PORT, new Prop()));
            add(new Property(Property.STRING_TYPE, "serviceRegion"));
            add(new Property(Property.STRING_TYPE, DbCredentials.CONNECTION_STRING,
                    DbCredentials.CONNECTION_STRING_DESCRIPTION, new Prop("textarea")));
            add(new Property(Property.BOOL_TYPE, DbCredentials.CACHE_SCHEMA));
            add(new Property(Property.BOOL_TYPE, DbCredentials.CACHE_RESULTS));
            add(new Property(Property.STRING_TYPE, DbCredentials.CACHE_INVALIDATE_SCHEDULE));
        }};
        descriptor.credentialsTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, "accessKey"));
            add(new Property(Property.STRING_TYPE, "secretAccessKey", new Prop("password")));
        }};
    }

    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = new java.util.Properties();
        if (!conn.hasCustomConnectionString() && conn.ssl()) {
            properties.setProperty("authScheme", "IAMSigV4");
            properties.setProperty("enableSsl", "true");
        }

        if (conn.credentials != null) {
            if (conn.credentials.parameters.get("accessKey") != null)
                System.setProperty("aws.accessKeyId",
                        conn.credentials.parameters.get("accessKey").toString());
            if (conn.credentials.parameters.get("secretAccessKey") != null)
                System.setProperty("aws.secretKey",
                        conn.credentials.parameters.get("secretAccessKey").toString());
            if (conn.credentials.parameters.get("serviceRegion") != null)
                properties.setProperty("serviceRegion",  conn.credentials.parameters.get("serviceRegion").toString());
        }

        return properties;
    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        prepareProvider();
        return CustomDriverManager.getConnection(getConnectionString(conn), getProperties(conn), driverClassName);
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ";port=" + conn.getPort();
        return "jdbc:neptune:gremlin://" + conn.getServer() + port;
    }

    @Override
    public boolean autoInterpolation() {
        return false;
    }

    @Override
    protected String interpolateBool(FuncParam param) {
        return ((boolean) param.value) ? "true" : "false";
    }
}

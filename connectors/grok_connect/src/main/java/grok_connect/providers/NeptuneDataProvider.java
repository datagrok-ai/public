package grok_connect.providers;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncParam;
import grok_connect.utils.*;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.net.URL;
import java.net.URLClassLoader;
import java.sql.Connection;
import java.sql.SQLException;
import java.util.*;

public class NeptuneDataProvider extends JdbcDataProvider {
    private static final List<String> SUPPORTED_VERSIONS =
            Collections.unmodifiableList(Arrays.asList("< 1.2.0.0 and >= 1.1.1.0", ">= 1.2.0.0"));
    public NeptuneDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "software.aws.neptune.NeptuneDriver";
        descriptor = new DataSource();
        descriptor.type = "Neptune";
        descriptor.commentStart = "//";
        descriptor.description = "Amazon Neptune is a fully managed graph database service";
        Property engine = new Property(Property.STRING_TYPE, "engineVersion");
        engine.choices = SUPPORTED_VERSIONS;
        descriptor.connectionTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.SERVER));
            add(new Property(Property.INT_TYPE, DbCredentials.PORT, new Prop()));
            add(new Property(Property.STRING_TYPE, "serviceRegion"));
            add(new Property(Property.STRING_TYPE, DbCredentials.CONNECTION_STRING,
                    DbCredentials.CONNECTION_STRING_DESCRIPTION, new Prop("textarea")));
            add(engine);
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
        }
        if (conn.parameters.get("serviceRegion") != null)
            properties.setProperty("serviceRegion",  conn.parameters.get("serviceRegion").toString());

        return properties;
    }

    public Connection getConnection(DataConnection conn) throws SQLException {
        String engineVersion = conn.get("engineVersion");
        if (engineVersion.equals(SUPPORTED_VERSIONS.get(0))) {
            loadDriver("neptune-jdbc-2.0.0-all.jar");
        } else {
            loadDriver("amazon-neptune-jdbc-driver-3.0.0.jar");
        }
        return CustomDriverManager.getConnection(getConnectionString(conn), getProperties(conn), driverClassName);
    }

    private void loadDriver(String name) {
       String baseDir = getBaseDirName();
        try {
            URLClassLoader child = new URLClassLoader (new URL[] {new URL(String.format("jar:file://%s/lib/%s!/",
                    baseDir, name)).toURI().toURL()},
                    NeptuneDataProvider.class.getClassLoader());
            Class.forName(driverClassName, true, child);
        } catch (MalformedURLException | URISyntaxException | ClassNotFoundException e) {
            throw new RuntimeException(e);
        }
    }

    private String getBaseDirName() {
        Properties props = GrokConnect.getInfo();
        return props.get("basedir").toString();
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

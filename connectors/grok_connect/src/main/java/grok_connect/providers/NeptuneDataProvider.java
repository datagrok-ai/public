package grok_connect.providers;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncParam;
import grok_connect.utils.*;
import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLClassLoader;
import java.sql.Connection;
import java.sql.Driver;
import java.sql.SQLException;
import java.util.*;

public class NeptuneDataProvider extends JdbcDataProvider {
    private static final Map<String, String> SUPPORTED_VERSIONS = new HashMap<>();
    private final Map<String, Driver> drivers;

    static {
        SUPPORTED_VERSIONS.put("< 1.2.0.0", "neptune-jdbc-2.0.0-all.jar");
        SUPPORTED_VERSIONS.put(">= 1.2.0.0", "amazon-neptune-jdbc-driver-3.0.0-all.jar");
    }

    public NeptuneDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "software.aws.neptune.NeptuneDriver";
        descriptor = new DataSource();
        descriptor.type = "Neptune";
        descriptor.commentStart = "//";
        descriptor.description = "Amazon Neptune is a fully managed graph database service";
        Property engine = new Property(Property.STRING_TYPE, DbCredentials.ENGINE_VERSION,
                "Amazon Neptune engine version");
        engine.choices = new ArrayList<>(SUPPORTED_VERSIONS.keySet());
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

        Driver driver1 = getDriver("neptune-jdbc-2.0.0-all.jar");
        Driver driver2 = getDriver("amazon-neptune-jdbc-driver-3.0.0-all.jar");
        drivers = new HashMap<>();
        drivers.put("< 1.2.0.0", driver1);
        drivers.put(">= 1.2.0.0", driver2);
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
        String engineVersion = conn.get(DbCredentials.ENGINE_VERSION);
        if (engineVersion == null) {
            throw new RuntimeException("Engine version is mandatory");
        }
        return drivers.get(engineVersion).connect(getConnectionString(conn), getProperties(conn));
    }

    private Driver getDriver(String name) {
        String baseDir = getBaseDirName();
        File file = new File( String.format("%s/lib/%s", baseDir, name));
        try {
            URLClassLoader child = new URLClassLoader (new URL[] {file.toURI().toURL()},
                    this.getClass().getClassLoader());
            Class<?> aClass = Class.forName(driverClassName, true, child);
            return (Driver) aClass.getConstructor().newInstance();
        } catch (MalformedURLException | ClassNotFoundException | InvocationTargetException
                 | InstantiationException |
                 IllegalAccessException | NoSuchMethodException e) {
            throw new RuntimeException("Something went wrong when getting driver", e);
        }
    }

    private String getBaseDirName() {
        Properties props = GrokConnect.properties;
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

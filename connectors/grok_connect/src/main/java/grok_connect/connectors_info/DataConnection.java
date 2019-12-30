package grok_connect.connectors_info;

import java.util.*;


public class DataConnection
{
    private static final String SERVER = "server";
    private static final String DB = "db";
    private static final String PORT = "port";

    public String id;
    public String dataSource;
    public String connectionString;
    public Credentials credentials;

    public String getServer() { return (String)parameters.get(SERVER); }
    public String getDb() { return (String)parameters.get(DB); }
    public String getPort() {
        Object port = parameters.get(PORT);
        return port == null ? null : String.valueOf(((Double)parameters.get(PORT)).intValue());
    }

    public Map<String, Object> parameters = new TreeMap<>(String.CASE_INSENSITIVE_ORDER);
    public Map<String, String> tags = new TreeMap<>(String.CASE_INSENSITIVE_ORDER);

    public String get(String parameter) {
        return (String)parameters.get(parameter);
    }

    public DataConnection() { }

    public DataConnection(String connectionString)
    {
        this.connectionString = connectionString;
    }
}

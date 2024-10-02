package grok_connect.connectors_info;

import java.util.*;
import grok_connect.utils.*;

public class DbCredentials
{
    public static final String ACCESS_SERVER = "accessServer";
    public static final String SERVER = "server";
    public static final String DOMAIN = "domain";
    public static final String KEYSPACE = "keySpace";
    public static final String DB = "db";
    public static final String INITIAL_CATALOG = "initialCatalog";
    public static final String SCHEMA = "schema";
    public static final String LOGIN = "login";
    public static final String PASSWORD = "password";

    public static final String PORT = "port";
    public static final String CONNECTION_STRING = "connString";
    public static final String SSL = "ssl";
    public static final String CACHE_SCHEMA = "cacheSchema";
    public static final String CACHE_RESULTS = "cacheResults";
    public static final String CACHE_INVALIDATE_SCHEDULE = "cacheInvalidateSchedule";
    public static final String DB_DESCRIPTION = "Database name";
    public static final String CONNECTION_STRING_DESCRIPTION = "When specified, this connection string overrides " +
            "all other parameters except 'login' and 'password'";
    public static final String ACCOUNT_LOCATOR = "accountLocator";
    public static final String REGION_ID = "region";
    public static final String CLOUD = "cloud";
    public static final String WAREHOUSE = "warehouse";
    public static final String ACCOUNT = "account";
    public static final String S3OutputLocation = "S3OutputLocation";
    public static final String S3OutputEncOption = "S3OutputEncOption";
    public static final String ACCESS_KEY = "accessKey";
    public static final String SECRET_KEY = "secretKey";
    public static final String VPC_ENDPOINT = "VPCEndpoint";
    public static final String ROLE = "role";
    public static final String UID = "UID";
    public static final String PWD = "PWD";

    public String server;
    public String port;
    public String db;
    public String login;
    public String password;

    public DbCredentials() { }

    public DbCredentials(String server, String port, String db, String login, String password)
    {
        this.server = server;
        this.port = port;
        this.db = db;
        this.login = login;
        this.password = password;
    }

    public static List<Property> dbConnectionTemplate = new ArrayList<Property>() {{
        add(new Property(Property.STRING_TYPE, DbCredentials.SERVER));
        add(new Property(Property.INT_TYPE, DbCredentials.PORT, new Prop()));
        add(new Property(Property.STRING_TYPE, DbCredentials.DB, DbCredentials.DB_DESCRIPTION));
        add(new Property(Property.STRING_TYPE, DbCredentials.CONNECTION_STRING,
                DbCredentials.CONNECTION_STRING_DESCRIPTION, new Prop("textarea")));
    }};

    public static List<Property> dbCredentialsTemplate = new ArrayList<Property>() {{
        add(new Property(Property.STRING_TYPE, DbCredentials.LOGIN));
        add(new Property(Property.STRING_TYPE, DbCredentials.PASSWORD, new Prop("password")));
    }};

    public DbCredentials(Map<String, Object> map)
    {
        server = (String)map.get(SERVER);
        db = (String)map.get(DB);
        login = (String)map.get(LOGIN);
        password = (String)map.get(PASSWORD);
    }
}

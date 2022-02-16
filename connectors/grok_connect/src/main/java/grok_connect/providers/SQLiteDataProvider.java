package grok_connect.providers;

import java.sql.*;
import java.util.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class SQLiteDataProvider extends JdbcDataProvider {
    public boolean autoInterpolation() {
        return false;
    }

    public SQLiteDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "org.sqlite.JDBC";

        descriptor = new DataSource();
        descriptor.type = "SQLite";
        descriptor.description = "Query SQLite database";
        descriptor.connectionTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.DB, DbCredentials.DB_DESCRIPTION));
            add(new Property(Property.STRING_TYPE, DbCredentials.CONNECTION_STRING,
                    DbCredentials.CONNECTION_STRING_DESCRIPTION, new Prop("textarea")));
        }};
        descriptor.credentialsTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.LOGIN));
            add(new Property(Property.STRING_TYPE, DbCredentials.PASSWORD, new Prop("password")));
        }};
    }

    public String getConnectionStringImpl(DataConnection conn) {
        return "jdbc:sqlite:" + conn.getDb();
    }
}

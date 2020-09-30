package grok_connect.providers;

import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.table_query.AggrFunctionInfo;
import grok_connect.table_query.Stats;
import grok_connect.utils.Property;
import serialization.Types;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class SnowflakeDataProvider extends JdbcDataProvider{

    public SnowflakeDataProvider() {
        descriptor = new DataSource();
        descriptor.type = "Snowflake";
        descriptor.description = "Query Snowflake database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.nameBrackets = "\"";

        //TODO: .*
        descriptor.typesMap = new HashMap<String, String>() {{
            put("number", Types.FLOAT);
            put("decimal", Types.FLOAT);
            put("numeric", Types.FLOAT);

            put("int", Types.INT);
            put("integer", Types.INT);
            put("bigint", Types.BIG_INT);
            put("smallint", Types.INT);

            put("float", Types.FLOAT);
            put("float4", Types.FLOAT);
            put("float8", Types.FLOAT);
            put("double", Types.FLOAT);
            put("double precision", Types.FLOAT);
            put("real", Types.FLOAT);

            put("varchar", Types.STRING);
            put("char", Types.STRING);
            put("character", Types.STRING);
            put("string", Types.STRING);
            put("text", Types.STRING);

            put("boolean", Types.BOOL);

            put("date", Types.DATE_TIME);
            put("datetime", Types.DATE_TIME);
            put("time", Types.DATE_TIME);
            put("timestamp", Types.DATE_TIME);

        }};
        descriptor.aggregations.add(new AggrFunctionInfo(Stats.STDEV, "stddev(#)", Types.dataFrameNumericTypes));

    }

    public Connection getConnection(DataConnection conn) throws ClassNotFoundException, SQLException {
        Class.forName("com.snowflake.client.jdbc.SnowflakeDriver");
        java.util.Properties properties = defaultConnectionProperties(conn);

        if (!conn.hasCustomConnectionString() && conn.ssl())
            properties.setProperty("ssl", "on");

        properties.put("db", conn.getDb());
        return DriverManager.getConnection(getConnectionString(conn), properties);
    }

    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:snowflake://" + conn.getServer() + port;
    }
}

package grok_connect.providers;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Proxy;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;

import grok_connect.table_query.AggrFunctionInfo;
import grok_connect.table_query.GroupAggregation;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;
import org.json.JSONObject;
import serialization.Column;
import serialization.DataFrame;
import serialization.StringColumn;
import serialization.Types;

import java.util.ArrayList;
import java.util.Base64;
import java.util.HashMap;

public class BigQueryDataProvider extends JdbcDataProvider {
    private static final String SERVICE_ACCOUNT_METHOD = "Service Account/OAuthType=0";
    private static final String OAUTH_METHOD = "OAuthType=2";

    public BigQueryDataProvider() {
        driverClassName = "com.simba.googlebigquery.jdbc42.Driver";

        descriptor = new DataSource();
        descriptor.type = "BigQuery";
        descriptor.description = "Query BigQuery database";
        descriptor.connectionTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.CONNECTION_STRING,
                    DbCredentials.CONNECTION_STRING_DESCRIPTION, new Prop("textarea")));
            add(new Property(Property.STRING_TYPE, "projectId", "ID of project"));
        }};
        descriptor.credentialsTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.OAUTH_SERVICE_ACCOUNT_EMAIL, "Set the OAuthServiceAcctEmail property to your Google service account email address.", SERVICE_ACCOUNT_METHOD));
            add(new Property(Property.STRING_TYPE, DbCredentials.SECRET_KEY, "Key file that is used to authenticate the\n" +
                    "service account email address. This parameter supports keys in .json format.", SERVICE_ACCOUNT_METHOD, new Prop("rsa"), ".json"));
            add(new Property(Property.STRING_TYPE, "#OAuthAccessToken", null, OAUTH_METHOD, new Prop("password")));
        }};

        descriptor.nameBrackets = "`";
        descriptor.typesMap = new HashMap<String, String>() {{
            put("string", Types.STRING);
            put("#int.*", Types.INT);
            put("#float.*", Types.FLOAT);
            put("numeric", Types.FLOAT);
            put("bignumeric", Types.BIG_INT);
            put("boolean", Types.BOOL);
            put("date", Types.DATE_TIME);
            put("datetime", Types.DATE_TIME);
            put("time", Types.DATE_TIME);
            put("#timestamp.*", Types.DATE_TIME);
            put("#interval.*", Types.STRING);
            put("array", Types.LIST);
            put("bytes", Types.BLOB);
        }};
        descriptor.canBrowseSchema = true;
    }

    @Override
    public void configureAutoCommit(Connection connection) {
    }

    @Override
    public Connection getConnection(DataConnection conn) throws SQLException, GrokConnectException {
        prepareProvider();
        if (conn.credentials == null)
            throw new GrokConnectException("Credentials can't be null");

        boolean resolvedByDatagrok = conn.credentials.parameters.getOrDefault("#dg_resolved", Boolean.FALSE).equals(true);
        com.simba.googlebigquery.jdbc42.DataSource ds = getDataSource(conn, resolvedByDatagrok);
        if (resolvedByDatagrok)
            return ds.getConnection();
        String key = getDataSourceKey(conn, ds);
        Connection connection = ConnectionPool.getConnection(ds, key);
        return (Connection) Proxy.newProxyInstance(
                conn.getClass().getClassLoader(),
                new Class[]{Connection.class},
                (proxy, method, args) -> {
                    try {
                        return method.invoke(connection, args);
                    } catch (InvocationTargetException e) {
                        Throwable target = e.getTargetException();
                        if (target instanceof SQLException)
                            throw wrapBigQuerySQLException((SQLException) target);
                        else
                            throw target;
                    }
                }
        );
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        return "jdbc:bigquery://https://www.googleapis.com/bigquery/v2:443;";
    }

    @Override
    public DataFrame getSchemas(DataConnection connection) throws GrokConnectException {
        try (Connection dbConnection = getConnection(connection);
             ResultSet schemas = dbConnection.getMetaData().getSchemas()) {
            DataFrame result = new DataFrame();
            Column<?> tableSchemaColumn = new StringColumn("table_schema");
            result.addColumn(tableSchemaColumn);
            while (schemas.next())
                result.addRow(schemas.getString(1));
            return result;
        } catch (SQLException e) {
            throw new GrokConnectException(e);
        }
    }

    @Override
    public DataFrame getSchema(DataConnection connection, String schema, String table, boolean includeKeyInfo) throws
            QueryCancelledByUser, GrokConnectException {
        try (Connection dbConnection = getConnection(connection);
             ResultSet columns = dbConnection.getMetaData().getColumns(null, schema, table, null)) {
            DataFrame result = DataFrame.fromColumns(new StringColumn("table_schema"),
                    new StringColumn("table_name"), new StringColumn("column_name"),
                    new StringColumn("data_type"));
            while (columns.next())
                result.addRow(columns.getString(2), columns.getString(3),
                        columns.getString(4), columns.getString(6));
            return result;
        } catch (SQLException e) {
            throw new GrokConnectException(e);
        }
    }

    @Override
    public String aggrToSql(GroupAggregation aggr) {
        AggrFunctionInfo funcInfo = null;
        for (AggrFunctionInfo info: descriptor.aggregations) {
            if (info.functionName.equals(aggr.aggType)) {
                funcInfo = info;
                break;
            }
        }
        if (funcInfo != null) {
            String brackets = descriptor.nameBrackets;
            String sql = GrokConnectUtil.isNotEmpty(aggr.colName) ? funcInfo.dbFunctionName.replaceAll("#", addBrackets(aggr.colName)) : funcInfo.dbFunctionName;
            return sql + " as " + (GrokConnectUtil.isNotEmpty(aggr.resultColName)
                    ?  brackets.charAt(0) + aggr.resultColName + brackets.charAt(brackets.length() - 1)
                    : brackets.charAt(0) + funcInfo.functionName + "_" + aggr.colName + brackets.charAt(brackets.length() - 1));
        }
        else
            return null;
    }

    private com.simba.googlebigquery.jdbc42.DataSource getDataSource(DataConnection conn, boolean resolvedByDatagrok) throws GrokConnectException {
        try {
            com.simba.googlebigquery.jdbc42.DataSource ds = new com.simba.googlebigquery.jdbc42.DataSource();
            ds.setLogDirectory("");
            ds.setURL(getConnectionString(conn));
            if (!conn.hasCustomConnectionString()) {
                ds.setProjectId(conn.get("projectId"));
                ds.setLogLevel("2");
            }

            String method = (String) conn.credentials.parameters.get("#chosen-auth-method");
            if (GrokConnectUtil.isEmpty(method) && !resolvedByDatagrok)
                throw new GrokConnectException("Authentication method was not set");
            if (resolvedByDatagrok || method.equals(OAUTH_METHOD)) {
                String token = (String) conn.credentials.parameters.get(DbCredentials.SECRET_KEY);
                if (GrokConnectUtil.isEmpty(token))
                    throw new GrokConnectException("Invalid OAuth token");
                ds.setOAuthType(2);
                ds.setOAuthAccessToken(token);
            }
            else if (method.equals(SERVICE_ACCOUNT_METHOD)) {
                String serviceEmail = (String) conn.credentials.parameters.get(DbCredentials.OAUTH_SERVICE_ACCOUNT_EMAIL);
                if (GrokConnectUtil.isEmpty(serviceEmail))
                    throw new GrokConnectException("OAuthServiceAcctEmail is mandatory for OAuthType=0");
                String jsonSecretKey = (String) conn.credentials.parameters.get(DbCredentials.SECRET_KEY);
                if (GrokConnectUtil.isEmpty(jsonSecretKey))
                    throw new GrokConnectException("Secret key file should be provided for OAuthType=0");
                ds.setOAuthType(0);
                ds.setOAuthServiceAcctEmail(serviceEmail);
                ds.setOAuthPvtKey(new String(Base64.getDecoder().decode(jsonSecretKey)));
            }
            else
                throw new GrokConnectException("Unsupported authentication method for BigQuery");

            return ds;
        } catch (Exception e) {
            throw new GrokConnectException(e);
        }
    }

    private String getDataSourceKey(DataConnection conn, com.simba.googlebigquery.jdbc42.DataSource ds) throws GrokConnectException {
        try {
            String method = (String) conn.credentials.parameters.get("#chosen-auth-method");
            String privateKeyHash;
            StringBuilder rawKey = new StringBuilder();
            rawKey.append("url=");
            rawKey.append(getConnectionString(conn));
            if (!conn.hasCustomConnectionString() && GrokConnectUtil.isNotEmpty(ds.getProjectId())) {
                rawKey.append("projectId=");
                rawKey.append(ds.getProjectId());
            }
            if (method.equals(SERVICE_ACCOUNT_METHOD)) {
                String serviceEmail = (String) conn.credentials.parameters.get(DbCredentials.OAUTH_SERVICE_ACCOUNT_EMAIL);
                rawKey.append("serviceEmail=");
                rawKey.append(serviceEmail);
                String jsonSecretKey = (String) conn.credentials.parameters.get(DbCredentials.SECRET_KEY);
                privateKeyHash = getKeyHash(jsonSecretKey);
            }
            else
                privateKeyHash = getKeyHash((String) conn.credentials.parameters.get(DbCredentials.SECRET_KEY));

            rawKey.append("key=");
            rawKey.append(privateKeyHash);
            return Base64.getEncoder().encodeToString(rawKey.toString().getBytes(StandardCharsets.UTF_8));
        } catch (NoSuchAlgorithmException e) {
            throw new GrokConnectException(e);
        }
    }

    private String getKeyHash(String key) throws NoSuchAlgorithmException {
        MessageDigest digest = MessageDigest.getInstance("SHA-256");
        return Base64.getEncoder().encodeToString(digest.digest(key.trim().getBytes(StandardCharsets.UTF_8)));
    }

    private SQLException wrapBigQuerySQLException(SQLException e) {
        String message = e.getMessage();
        String newMessage = message;

        int jsonStart = message.indexOf("{");
        if (jsonStart != -1) {
            String jsonPart = message.substring(jsonStart);
            try {
                JSONObject obj = new JSONObject(jsonPart);
                if (obj.has("message")) {
                    newMessage = obj.getString("message");
                }
            } catch (Exception ex) {
                // Parsing failed, keep original message
            }
        }

        SQLException wrapped = new SQLException(newMessage, e.getSQLState(), e.getErrorCode(), e);
        wrapped.setStackTrace(e.getStackTrace());
        return wrapped;
    }

    @Override
    public String getCommentsQuery(DataConnection connection) throws GrokConnectException {
        Object projectId = connection.get("projectId");
        if (projectId == null)
            throw new GrokConnectException("Project id should be provided to get comments");

        return String.format(
                "--input: string schema\n" +
                        "-- TABLE and VIEW descriptions\n" +
                        "(\n" +
                        "  SELECT\n" +
                        "      table_schema AS table_schema,\n" +
                        "      table_name AS object_name,\n" +
                        "      table_type AS object_type,\n" +
                        "      NULL AS column_name,\n" +
                        "      table_catalog,\n" +
                        "      table_schema,\n" +
                        "      table_name,\n" +
                        "      ddl AS full_ddl,           -- sometimes contains COMMENT\n" +
                        "      option_value AS comment\n" +
                        "  FROM `%s.@schema.INFORMATION_SCHEMA.TABLE_OPTIONS`\n" +
                        "  WHERE option_name = 'description'\n" +
                        ")\n" +
                        "\n" +
                        "UNION ALL\n" +
                        "\n" +
                        "-- COLUMN descriptions\n" +
                        "(\n" +
                        "  SELECT\n" +
                        "      table_schema AS table_schema,\n" +
                        "      table_name AS object_name,\n" +
                        "      'COLUMN' AS object_type,\n" +
                        "      column_name AS column_name,\n" +
                        "      NULL AS table_catalog,\n" +
                        "      table_schema,\n" +
                        "      table_name,\n" +
                        "      NULL AS full_ddl,\n" +
                        "      description AS comment\n" +
                        "  FROM `%s.@schema.INFORMATION_SCHEMA.COLUMN_FIELD_PATHS`\n" +
                        ")\n" +
                        "\n" +
                        "ORDER BY object_type, object_name, column_name;\n",
                projectId,
                projectId
        );
    }
}

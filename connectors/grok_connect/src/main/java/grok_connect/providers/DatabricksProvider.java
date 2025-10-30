package grok_connect.providers;

import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.utils.*;
import serialization.DataFrame;
import serialization.IntColumn;
import serialization.StringColumn;
import serialization.Types;
import java.sql.Connection;
import java.sql.DatabaseMetaData;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.*;

public class DatabricksProvider extends JdbcDataProvider {
    private static final String PAT_METHOD = "Personal Access Token";
    private static final String OAUTH_METHOD = "OAuth Client Credentials";

    public DatabricksProvider() {
        driverClassName = "com.databricks.client.jdbc.Driver";

        descriptor = new DataSource();
        descriptor.type = "Databricks";
        descriptor.description = "Query Databricks warehouse";
        descriptor.requiresFullyQualifiedTable = true;
        descriptor.connectionTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.CONNECTION_STRING,
                    DbCredentials.CONNECTION_STRING_DESCRIPTION, new Prop("textarea")));
            add(new Property(Property.STRING_TYPE, "workspaceURL", "Your Databricks workspace host. You can copy this from your browser’s address bar — for example dbc-52240d8b-a70a.cloud.databricks.com."));
            add(new Property(Property.STRING_TYPE, "httpPath", "The unique path of your SQL Warehouse or endpoint. You can find it in Databricks → Compute → SQL Warehouses → Connection details."));
            add(new Property(Property.STRING_TYPE, "Catalog", "Optional. Unity Catalog name that contains your data (for example samples or main). If omitted, Databricks defaults to main."));
        }};

        descriptor.credentialsTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, "token",
                    "Your Databricks personal access token (PAT). Generate one in Settings → Developer → Access tokens. It acts as your password.", PAT_METHOD, new Prop("password")));
            add(new Property(Property.STRING_TYPE, "clientId",
                    "The application (client) ID of your Azure AD or Databricks service principal.", OAUTH_METHOD));
            add(new Property(Property.STRING_TYPE, "clientSecret",
                    "The secret associated with the client ID.", OAUTH_METHOD, new Prop("password")));
            add(new Property(Property.STRING_TYPE, "tenantId",
                    "Azure Directory ID (tenant GUID) used for OAuth authentication. Leave empty for AWS or GCP workspaces.", OAUTH_METHOD));
        }};

        descriptor.nameBrackets = "`";
        descriptor.defaultSchema = "default";
        descriptor.canBrowseSchema = true;
        descriptor.typesMap = new HashMap<String, String>() {{
            put("tinyint", Types.INT);
            put("smallint", Types.INT);
            put("integer", Types.INT);
            put("int", Types.INT);
            put("bigint", Types.BIG_INT);
            put("real", Types.FLOAT);
            put("double", Types.FLOAT);
            put("#decimal.*", Types.FLOAT);
            put("boolean", Types.BOOL);
            put("#char.*", Types.STRING);
            put("#varchar.*", Types.STRING);
            put("string", Types.STRING);
            put("text", Types.STRING);
            put("date", Types.DATE_TIME);
            put("#timestamp.*", Types.DATE_TIME);
            put("#varbinary.*", Types.BLOB);
        }};
    }

    @Override
    public Properties getProperties(DataConnection conn) {
        Properties properties = new Properties();
        if (!conn.hasCustomConnectionString()) {
            properties.setProperty("SSL", "1"); // Always true
            setIfNotEmpty(properties, "ConnCatalog", (String) conn.get("Catalog"));
        }

        String method = (String) conn.credentials.parameters.get("#chosen-auth-method");

        if (GrokConnectUtil.isNotEmpty(method) && method.equals(OAUTH_METHOD)) {
            setIfNotNull(properties, DbCredentials.OAUTH2_CLIENT_ID, (String) conn.credentials.parameters.get("clientId"));
            setIfNotNull(properties, DbCredentials.OAUTH2_SECRET, (String) conn.credentials.parameters.get("clientSecret"));
            setIfNotEmpty(properties, "AzureTenantId", (String) conn.credentials.parameters.get("tenantId"));
            properties.setProperty("AuthMech", "11");
            properties.setProperty("Auth_Flow", "1");
        }
        else {
            properties.setProperty("AuthMech", "3");
            properties.setProperty(DbCredentials.UID, "token");
            setIfNotNull(properties, DbCredentials.PWD, (String) conn.credentials.parameters.get("token"));
        }

        return properties;
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        String workspaceUrl = conn.get("workspaceURL");
        if (GrokConnectUtil.isEmpty(workspaceUrl)) {
            workspaceUrl = workspaceUrl.replaceFirst("(?i)^https?://", "");
            int slashIndex = workspaceUrl.indexOf('/');
            if (slashIndex != -1)
                workspaceUrl = workspaceUrl.substring(0, slashIndex);
            int colonIndex = workspaceUrl.indexOf(':');
            if (colonIndex != -1)
                workspaceUrl = workspaceUrl.substring(0, colonIndex);
        }
        String url = "jdbc:databricks://" + workspaceUrl + ":443" +  "/" + "default";
        String httpPath = conn.get("httpPath");
        if (GrokConnectUtil.isNotEmpty(httpPath))
            url += ";httpPath=" + httpPath + ";";
        return url;
    }

    @Override
    public DataFrame getSchemas(DataConnection connection) throws QueryCancelledByUser, GrokConnectException {
        try (Connection dbConnection = getConnection(connection);
                 ResultSet schemas = dbConnection.getMetaData().getSchemas()) {
            DataFrame result = DataFrame.fromColumns(new StringColumn("table_schema"));
            String db = connection.get("Catalog");
            while (schemas.next()) {
                String catalog = schemas.getString(2);
                if (GrokConnectUtil.isEmpty(db) || catalog == null || catalog.equals(db))
                    result.addRow(schemas.getString(1));
            }
            return result;
        } catch (SQLException e) {
            throw new GrokConnectException(e);
        }
    }

    @Override
    public DataFrame getSchema(DataConnection connection, String schema, String table) throws
            QueryCancelledByUser, GrokConnectException {

        try (Connection dbConnection = getConnection(connection)) {
            DataFrame result = DataFrame.fromColumns(new StringColumn("table_catalog"), new StringColumn("table_schema"),
                    new StringColumn("table_name"), new StringColumn("column_name"),
                    new StringColumn("data_type"), new IntColumn("is_view"));
            DatabaseMetaData metaData = dbConnection.getMetaData();
            String db = GrokConnectUtil.isEmpty(connection.get("Catalog")) ? null : connection.get("Catalog");
            try (ResultSet tables = metaData.getTables(db, schema, "%", new String[] {"TABLE", "VIEW"})) {
                while (tables.next()) {
                    String tableCatalog = tables.getString("TABLE_CAT");
                    String tableSchema = tables.getString("TABLE_SCHEM");
                    String tableName = tables.getString("TABLE_NAME");
                    Integer isView = "VIEW".equals(tables.getString("TABLE_TYPE")) ? 1 : 0;
                    boolean hasColumns = false;
                    try (ResultSet columns = metaData.getColumns(tableCatalog, tableSchema, tableName, "%")) {
                        while (columns.next()) {
                            hasColumns = true;
                            String columnName = columns.getString("COLUMN_NAME");
                            String dataType = columns.getString("TYPE_NAME");
                            result.addRow(tableCatalog, tableSchema, tableName, columnName, dataType, isView);
                        }
                    }
                    if (!hasColumns)
                        result.addRow(tableCatalog, tableSchema, tableName, null, null, isView);
                }
            }
            return result;
        } catch (SQLException e) {
            throw new GrokConnectException(e);
        }
    }

    @Override
    public DataFrame getForeignKeys(DataConnection conn, String schema) {
        // Return empty df since foreign keys are not supported
        return DataFrame.fromColumns(new StringColumn("table_schema"),
                new StringColumn("constraint_name"), new StringColumn("table_name"),
                new StringColumn("column_name"), new StringColumn("foreign_table_name"), new StringColumn("foreign_column_name"));
    }

    public static String simplifyDatabricksError(String msg) {
        if (GrokConnectUtil.isEmpty(msg))
            return msg;

        msg = msg.replace("\r\n", "\n");

        msg = msg.replaceAll("(?s)(\\n\\s*at org\\.apache\\.spark[^\n]*)[\\s\\S]*", "");
        msg = msg.replaceAll("(?s)(\\nCaused by: org\\.apache\\.spark[^\n]*)[\\s\\S]*", "");

        msg = msg.trim().replaceAll("\\n{2,}", "\n");

        return msg;
    }
}

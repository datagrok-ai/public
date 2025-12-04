package grok_connect.providers;

import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.utils.*;
import serialization.DataFrame;
import serialization.IntColumn;
import serialization.StringColumn;
import serialization.Types;

import java.sql.*;
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
            if (!conn.hasCustomConnectionString())
                setIfNotEmpty(properties, "ConnCatalog", conn.get("Catalog"));
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
    public DataFrame getSchemas(DataConnection connection)
            throws QueryCancelledByUser, GrokConnectException {

        try (Connection db = getConnection(connection)) {

            DataFrame result =
                    DataFrame.fromColumns(new StringColumn("table_schema"));

            String catalog = getCatalog(connection, db);

            String infoSchemaPath = getInfoSchemaPath(db, catalog);

            String sql =
                    "SELECT schema_name " +
                            "FROM " + infoSchemaPath + ".schemata " +
                            "ORDER BY schema_name";

            try (PreparedStatement ps = db.prepareStatement(sql);
                 ResultSet rs = ps.executeQuery()) {

                while (rs.next())
                    result.addRow(rs.getString(1));
            }

            return result;

        } catch (SQLException e) {
            throw new GrokConnectException(simplifyDatabricksError(e.getMessage()));
        }
    }

    private String getInfoSchemaPath(Connection db, String catalog) throws SQLException {
        return isUnityCatalog(db) ? addBrackets(catalog) + ".information_schema" : "information_schema";
    }

    @Override
    public DataFrame getSchema(DataConnection connection, String schema, String table, boolean includeKeyInfo)
            throws QueryCancelledByUser, GrokConnectException {

        try (Connection db = getConnection(connection)) {

            DataFrame result = DataFrame.fromColumns(
                    new StringColumn("table_catalog"),
                    new StringColumn("table_schema"),
                    new StringColumn("table_name"),
                    new StringColumn("column_name"),
                    new StringColumn("data_type"),
                    new IntColumn("is_view")
            );

            String catalog = getCatalog(connection, db);
            String infoSchemaPath = getInfoSchemaPath(db, catalog);
            StringBuilder sql = new StringBuilder(
                    "SELECT t.table_catalog, t.table_schema, t.table_name, " +
                            "       c.column_name, c.data_type, " +
                            "       CASE WHEN t.table_type = 'VIEW' THEN 1 ELSE 0 END AS is_view " +
                            "FROM " + infoSchemaPath + ".tables t " +
                            "LEFT JOIN " + infoSchemaPath + ".columns c " +
                            "  ON t.table_catalog = c.table_catalog " +
                            " AND t.table_schema  = c.table_schema " +
                            " AND t.table_name    = c.table_name "
            );


            sql.append("WHERE 1=1 ");

            if (GrokConnectUtil.isNotEmpty(schema))
                sql.append("AND LOWER(t.table_schema) = LOWER(?) ");

            if (GrokConnectUtil.isNotEmpty(table))
                sql.append("AND t.table_name = ? ");

            sql.append("ORDER BY t.table_schema, t.table_name, c.ordinal_position");

            try (PreparedStatement ps = db.prepareStatement(sql.toString())) {

                int idx = 1;
                if (schema != null)
                    ps.setString(idx++, schema);

                if (table != null)
                    ps.setString(idx, table);

                try (ResultSet rs = ps.executeQuery()) {
                    while (rs.next()) {
                        result.addRow(
                                rs.getString("table_catalog"),
                                rs.getString("table_schema"),
                                rs.getString("table_name"),
                                rs.getString("column_name"),
                                rs.getString("data_type"),
                                rs.getInt("is_view")
                        );
                    }
                }
            }

            return result;

        } catch (SQLException e) {
            throw new GrokConnectException(simplifyDatabricksError(e.getMessage()));
        }
    }

    private String getCatalog(DataConnection connection, Connection db) throws SQLException {
        String catalog = connection.hasCustomConnectionString() ? null : connection.get("Catalog");
        if (GrokConnectUtil.isEmpty(catalog)) {
            try (Statement st = db.createStatement();
                 ResultSet rs = st.executeQuery("SELECT current_catalog()")) {
                if (rs.next())
                    catalog = rs.getString(1);
                else
                    catalog = "hive_metastore";
            }
        }
        return catalog;
    }

    private boolean isUnityCatalog(Connection db) throws SQLException {
        // --- 1. Try SHOW CATALOGS ---
        try (Statement st = db.createStatement();
             ResultSet rs = st.executeQuery("SHOW CATALOGS")) {

            int count = 0;
            boolean hasNonHiveCatalog = false;

            while (rs.next()) {
                count++;
                String cat = rs.getString(1);
                if (!"hive_metastore".equalsIgnoreCase(cat))
                    hasNonHiveCatalog = true;
            }

            // Unity Catalog is enabled when:
            // - More than one catalog exists, OR
            // - At least one catalog is not hive_metastore
            if (count > 1 || hasNonHiveCatalog)
                return true;
        }

        // --- 2. Unity Catalog disabled always returns a pseudo "hive_metastore"
        // but this catalog does NOT expose a nested information_schema.
        // Test if <catalog>.information_schema works.
        try (Statement st = db.createStatement();
             ResultSet test = st.executeQuery(
                     "SELECT * FROM hive_metastore.information_schema.tables LIMIT 1")) {

            // If the above query works, this "hive_metastore" is a REAL UC catalog.
            return true;

        } catch (SQLException ignored) {
            // Not a real catalog → no UC
        }

        // --- 3. Additional fallback: UC info schema exists only under catalogs ---
        // Try probing system catalog which only exists under UC
        try (Statement st = db.createStatement();
             ResultSet rs = st.executeQuery("SHOW SCHEMAS IN system")) {

            // If this works: Unity Catalog
            return true;

        } catch (SQLException ignored) {
            // system catalog not found → UC disabled
        }

        // If none of the above heuristics indicate UC, treat as Hive-only
        return false;
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

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

@SuppressWarnings({"SqlDialectInspection", "SqlNoDataSourceInspection"})
public class DatabricksProvider extends JdbcDataProvider {
    private static final String PAT_METHOD = "Personal Access Token";
    private static final String OAUTH_METHOD = "OAuth Client Credentials";
    private static final String FEDERATED_SSO = "Federated (SSO)";

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
            add(new Property(Property.STRING_TYPE, DbCredentials.DB, "Optional. Unity Catalog name that contains your data (for example samples or main). If omitted, Databricks defaults to main.", null, "Catalog"));
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
            add(new Property(Property.STRING_TYPE, "#token", null, FEDERATED_SSO, new Prop("password")));
        }};

        descriptor.nameBrackets = "`";
        descriptor.defaultSchema = "default";
        descriptor.canBrowseSchema = true;
        descriptor.supportCatalogs = true;
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
                setIfNotEmpty(properties, "ConnCatalog", conn.getDb());
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
            String token = (String) (method.equals(FEDERATED_SSO) ? conn.credentials.parameters.get("#token") : conn.credentials.parameters.get("token"));
            properties.setProperty("AuthMech", "3");
            properties.setProperty(DbCredentials.UID, "token");
            setIfNotNull(properties, DbCredentials.PWD, token);
        }

        return properties;
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        String workspaceUrl = conn.get("workspaceURL");
        if (GrokConnectUtil.isNotEmpty(workspaceUrl)) {
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

            SQLException primaryError;
            try {
                return getSchemasFromInfoSchema(db, infoSchemaPath, catalog);
            } catch (SQLException e) {
                primaryError = e;
                logger.debug("Primary getSchemas query failed: {}", e.getMessage());
            }

            // Fallback: try the opposite information_schema path
            String fallbackPath;
            if (infoSchemaPath.contains("."))
                fallbackPath = "information_schema";
            else if (GrokConnectUtil.isNotEmpty(catalog))
                fallbackPath = addBrackets(catalog) + ".information_schema";
            else
                fallbackPath = null;

            if (fallbackPath != null) {
                try {
                    logger.debug("Fallback getSchemas query with path: {}", fallbackPath);
                    return getSchemasFromInfoSchema(db, fallbackPath, catalog);
                } catch (SQLException e) {
                    logger.debug("Fallback getSchemas query also failed: {}", e.getMessage());
                }
            }

            // Final fallback: use SHOW SCHEMAS which works on both systems
            if (GrokConnectUtil.isNotEmpty(catalog)) {
                try (Statement st = db.createStatement();
                     ResultSet rs = st.executeQuery("SHOW SCHEMAS IN " + addBrackets(catalog))) {

                    logger.debug("SHOW SCHEMAS fallback succeeded");
                    while (rs.next())
                        result.addRow(rs.getString(1));
                    return result;
                } catch (SQLException e) {
                    logger.debug("SHOW SCHEMAS fallback failed: {}", e.getMessage());
                }
            }

            // Last resort: SHOW SCHEMAS without catalog
            try (Statement st = db.createStatement();
                 ResultSet rs = st.executeQuery("SHOW SCHEMAS")) {

                logger.debug("SHOW SCHEMAS (no catalog) fallback succeeded");
                while (rs.next())
                    result.addRow(rs.getString(1));
                return result;
            } catch (SQLException e) {
                logger.debug("SHOW SCHEMAS (no catalog) fallback failed: {}", e.getMessage());
            }

            // All attempts failed, throw the original error
            throw primaryError;

        } catch (SQLException e) {
            throw new GrokConnectException(simplifyDatabricksError(e.getMessage()));
        }
    }

    private DataFrame getSchemasFromInfoSchema(Connection db, String infoSchemaPath, String catalog)
            throws SQLException {
        DataFrame result = DataFrame.fromColumns(new StringColumn("table_schema"));

        String sql;
        if (GrokConnectUtil.isNotEmpty(catalog)) {
            sql = "SELECT DISTINCT schema_name FROM " + infoSchemaPath + ".schemata " +
                    "WHERE catalog_name = ? ORDER BY schema_name";
        } else {
            sql = "SELECT DISTINCT schema_name FROM " + infoSchemaPath + ".schemata " +
                    "ORDER BY schema_name";
        }

        try (PreparedStatement ps = db.prepareStatement(sql)) {
            if (GrokConnectUtil.isNotEmpty(catalog))
                ps.setString(1, catalog);
            try (ResultSet rs = ps.executeQuery()) {
                while (rs.next())
                    result.addRow(rs.getString(1));
            }
        }
        return result;
    }

    private String getInfoSchemaPath(Connection db, String catalog) throws SQLException {
        Boolean isUnity = isUnityCatalog(db);
        if (isUnity == null) {
            // Detection inconclusive - try Unity Catalog path first if catalog available, will fallback if it fails
            if (GrokConnectUtil.isNotEmpty(catalog))
                return addBrackets(catalog) + ".information_schema";
            return "information_schema";
        }
        if (isUnity && GrokConnectUtil.isNotEmpty(catalog))
            return addBrackets(catalog) + ".information_schema";
        return "information_schema";
    }

    @Override
    public DataFrame getSchema(DataConnection connection, String schema, String table, boolean includeKeyInfo)
            throws QueryCancelledByUser, GrokConnectException {

        try (Connection db = getConnection(connection)) {

            String catalog = getCatalog(connection, db);
            String infoSchemaPath = getInfoSchemaPath(db, catalog);

            // Try primary path
            DataFrame result = tryGetSchemaWithPath(db, infoSchemaPath, catalog, schema, table);
            if (result != null)
                return result;

            // Fallback: try the opposite information_schema path
            String fallbackPath;
            if (infoSchemaPath.contains(".")) {
                // Was Unity path, try legacy
                fallbackPath = "information_schema";
            } else if (GrokConnectUtil.isNotEmpty(catalog)) {
                // Was legacy, try Unity (only if catalog is available)
                fallbackPath = addBrackets(catalog) + ".information_schema";
            } else {
                fallbackPath = null;
            }

            if (fallbackPath != null) {
                logger.debug("Primary getSchema query failed, trying fallback path: {}", fallbackPath);
                result = tryGetSchemaWithPath(db, fallbackPath, catalog, schema, table);
                if (result != null)
                    return result;
            }

            // Final fallback: use SHOW TABLES and DESCRIBE for each table
            logger.debug("Information schema queries failed, using SHOW TABLES fallback");
            return getSchemaViaShowTables(db, catalog, schema, table);

        } catch (SQLException e) {
            throw new GrokConnectException(simplifyDatabricksError(e.getMessage()));
        }
    }

    private DataFrame tryGetSchemaWithPath(Connection db, String infoSchemaPath, String catalog,
                                           String schema, String table) {
        // Try with table_catalog first
        DataFrame result = tryGetSchemaWithCatalogColumn(db, infoSchemaPath, catalog, schema, table, true);
        if (result != null)
            return result;

        // If failed, try without table_catalog (some legacy/driver configs don't have it)
        logger.debug("Retrying getSchema without table_catalog column");
        return tryGetSchemaWithCatalogColumn(db, infoSchemaPath, catalog, schema, table, false);
    }

    private DataFrame tryGetSchemaWithCatalogColumn(Connection db, String infoSchemaPath, String catalog,
                                                    String schema, String table, boolean useCatalogColumn) {
        DataFrame result = DataFrame.fromColumns(
                new StringColumn("table_catalog"),
                new StringColumn("table_schema"),
                new StringColumn("table_name"),
                new StringColumn("column_name"),
                new StringColumn("data_type"),
                new IntColumn("is_view")
        );

        StringBuilder sql = new StringBuilder();
        sql.append("SELECT ");
        if (useCatalogColumn)
            sql.append("t.table_catalog, ");
        sql.append("t.table_schema, t.table_name, c.column_name, c.data_type, ")
           .append("CASE WHEN t.table_type = 'VIEW' THEN 1 ELSE 0 END AS is_view ")
           .append("FROM ").append(infoSchemaPath).append(".tables t ")
           .append("LEFT JOIN ").append(infoSchemaPath).append(".columns c ")
           .append("  ON LOWER(t.table_schema) = LOWER(c.table_schema) ")
           .append(" AND LOWER(t.table_name) = LOWER(c.table_name) ");
        if (useCatalogColumn)
            sql.append(" AND LOWER(t.table_catalog) = LOWER(c.table_catalog) ");

        sql.append("WHERE 1=1 ");

        // Filter by catalog only if column exists and catalog is provided
        if (useCatalogColumn && GrokConnectUtil.isNotEmpty(catalog))
            sql.append("AND t.table_catalog = ? ");

        if (GrokConnectUtil.isNotEmpty(schema))
            sql.append("AND LOWER(t.table_schema) = LOWER(?) ");

        if (GrokConnectUtil.isNotEmpty(table))
            sql.append("AND t.table_name = ? ");

        sql.append("ORDER BY t.table_schema, t.table_name, c.ordinal_position");

        try (PreparedStatement ps = db.prepareStatement(sql.toString())) {
            int idx = 1;
            if (useCatalogColumn && GrokConnectUtil.isNotEmpty(catalog))
                ps.setString(idx++, catalog);

            if (GrokConnectUtil.isNotEmpty(schema))
                ps.setString(idx++, schema);

            if (GrokConnectUtil.isNotEmpty(table))
                ps.setString(idx, table);

            try (ResultSet rs = ps.executeQuery()) {
                while (rs.next()) {
                    result.addRow(
                            useCatalogColumn ? rs.getString("table_catalog") : catalog,
                            rs.getString("table_schema"),
                            rs.getString("table_name"),
                            rs.getString("column_name"),
                            rs.getString("data_type"),
                            rs.getInt("is_view")
                    );
                }
            }
            return result;
        } catch (SQLException e) {
            logger.debug("getSchema query failed with path {} (useCatalogColumn={}): {}",
                    infoSchemaPath, useCatalogColumn, e.getMessage());
            return null;
        }
    }

    private DataFrame getSchemaViaShowTables(Connection db, String catalog, String schema, String table)
            throws SQLException {
        DataFrame result = DataFrame.fromColumns(
                new StringColumn("table_catalog"),
                new StringColumn("table_schema"),
                new StringColumn("table_name"),
                new StringColumn("column_name"),
                new StringColumn("data_type"),
                new IntColumn("is_view")
        );

        List<String> schemas = new ArrayList<>();
        if (GrokConnectUtil.isNotEmpty(schema)) {
            schemas.add(schema);
        } else {
            String showSchemasSql = GrokConnectUtil.isNotEmpty(catalog)
                    ? "SHOW SCHEMAS IN " + addBrackets(catalog) : "SHOW SCHEMAS";
            try (Statement st = db.createStatement();
                 ResultSet rs = st.executeQuery(showSchemasSql)) {
                while (rs.next())
                    schemas.add(rs.getString(1));
            }
        }

        for (String targetSchema : schemas) {
            String showTablesSql = GrokConnectUtil.isNotEmpty(catalog)
                    ? "SHOW TABLES IN " + addBrackets(catalog) + "." + addBrackets(targetSchema)
                    : "SHOW TABLES IN " + addBrackets(targetSchema);

            List<String> tables = new ArrayList<>();
            try (Statement st = db.createStatement();
                 ResultSet rs = st.executeQuery(showTablesSql)) {
                while (rs.next()) {
                    String tableName = rs.getString("tableName");
                    if (table == null || table.equals(tableName))
                        tables.add(tableName);
                }
            }

            for (String tableName : tables) {
                String describeSql = GrokConnectUtil.isNotEmpty(catalog)
                        ? "DESCRIBE " + addBrackets(catalog) + "." + addBrackets(targetSchema) + "." + addBrackets(tableName)
                        : "DESCRIBE " + addBrackets(targetSchema) + "." + addBrackets(tableName);

                try (Statement st = db.createStatement();
                     ResultSet rs = st.executeQuery(describeSql)) {
                    while (rs.next()) {
                        String colName = rs.getString("col_name");
                        if (colName == null || colName.isEmpty() || colName.startsWith("#"))
                            continue;
                        result.addRow(catalog, targetSchema, tableName, colName, rs.getString("data_type"), 0);
                    }
                } catch (SQLException e) {
                    logger.debug("DESCRIBE failed for table {}: {}", tableName, e.getMessage());
                }
            }
        }

        return result;
    }

    private String getCatalog(DataConnection connection, Connection db) throws SQLException {
        String catalog = connection.getDb();
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

    /**
     * Detects if the connection uses Unity Catalog.
     * @return true if Unity Catalog, false if legacy Hive metastore, null if detection is inconclusive
     */
    private Boolean isUnityCatalog(Connection db) {
        // Method 1: Check SHOW CATALOGS - Unity Catalog has multiple catalogs or 'system' catalog
        try (Statement st = db.createStatement();
             ResultSet rs = st.executeQuery("SHOW CATALOGS")) {
            Set<String> catalogs = new HashSet<>();
            while (rs.next()) {
                catalogs.add(rs.getString(1).toLowerCase());
            }
            // Unity Catalog always has 'system' catalog; legacy Hive metastore does not
            if (catalogs.contains("system")) {
                logger.debug("Unity Catalog detected: 'system' catalog present");
                return true;
            }
            // If only hive_metastore exists and no system catalog, it's legacy
            if (catalogs.size() == 1 && catalogs.contains("hive_metastore")) {
                logger.debug("Legacy Hive metastore detected: only hive_metastore catalog present");
                return false;
            }
            // Multiple catalogs without 'system' is unusual but treat as Unity Catalog
            if (catalogs.size() > 1) {
                logger.debug("Unity Catalog detected: multiple catalogs present ({})", catalogs);
                return true;
            }
        } catch (SQLException e) {
            logger.debug("SHOW CATALOGS failed (may lack permissions): {}", e.getMessage());
        }

        // Method 2: Try to access system.information_schema (Unity Catalog specific)
        try (Statement st = db.createStatement();
             ResultSet ignored = st.executeQuery("SELECT 1 FROM system.information_schema.catalogs LIMIT 1")) {
            logger.debug("Unity Catalog detected via system.information_schema.catalogs");
            return true;
        } catch (SQLException e) {
            String msg = e.getMessage().toLowerCase();
            // If error indicates catalog/schema doesn't exist, it's likely legacy
            if (msg.contains("catalog") && (msg.contains("not found") || msg.contains("does not exist"))) {
                logger.debug("Legacy Hive metastore detected: system catalog does not exist");
                return false;
            }
            // Permission error - detection inconclusive
            logger.debug("system.information_schema query failed (may lack permissions): {}", e.getMessage());
        }

        // Method 3: Try legacy information_schema directly (no catalog prefix)
        // Note: This can work on both legacy AND Unity Catalog (with default catalog context),
        // so we can't definitively determine the system type - return null (inconclusive)
        try (Statement st = db.createStatement();
             ResultSet ignored = st.executeQuery("SELECT 1 FROM information_schema.tables LIMIT 1")) {
            logger.debug("Unqualified information_schema accessible - could be legacy or Unity with default catalog");
            return null;
        } catch (SQLException e) {
            logger.debug("Legacy information_schema query failed: {}", e.getMessage());
        }

        logger.debug("Unity Catalog detection inconclusive - will try Unity path with fallback");
        return null; // Inconclusive - caller should handle fallback
    }




    @Override
    public DataFrame getCatalogs(DataConnection connection) throws GrokConnectException {
        try (Connection db = getConnection(connection)) {
            DataFrame result = DataFrame.fromColumns(new StringColumn("catalog_name"));
            try (Statement st = db.createStatement();
                 ResultSet rs = st.executeQuery("SHOW CATALOGS")) {
                List<String> catalogs = new ArrayList<>();
                while (rs.next())
                    catalogs.add(rs.getString(1));
                Collections.sort(catalogs);
                for (String catalog : catalogs)
                    result.addRow(catalog);
                return result;
            } catch (SQLException e) {
                logger.debug("SHOW CATALOGS failed, falling back to JDBC metadata: {}", e.getMessage());
            }
            return readCatalogsFromMetadata(db);
        } catch (SQLException e) {
            throw new GrokConnectException(simplifyDatabricksError(e.getMessage()));
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

    @Override
    public String getCommentsQuery(DataConnection connection) throws GrokConnectException {
        try (Connection conn = getConnection(connection)) {
            Boolean isUnity = isUnityCatalog(conn);
            if (Boolean.FALSE.equals(isUnity))
                throw new GrokConnectException("Cannot get comments for legacy metastore");
            // If isUnity is true or null (inconclusive), proceed - query will fail naturally if not Unity
        } catch (SQLException e) {
            throw new GrokConnectException(e);
        }

        return "--input: string schema\n" +
                "-- SCHEMA (database.schema) comments\n" +
                "(\n" +
                "    SELECT \n" +
                "        s.catalog_name,\n" +
                "        s.schema_name as table_schema,\n" +
                "        'schema' AS object_type,\n" +
                "        NULL AS table_name,\n" +
                "        NULL AS column_name,\n" +
                "        s.comment AS comment\n" +
                "    FROM system.information_schema.schemata s\n" +
                "    WHERE s.schema_name = @schema\n" +
                ")\n" +
                "\n" +
                "UNION ALL\n" +
                "\n" +
                "-- TABLE comments\n" +
                "(\n" +
                "    SELECT \n" +
                "        t.table_catalog AS catalog_name,\n" +
                "        t.table_schema  AS table_schema,\n" +
                "        'table' AS object_type,\n" +
                "        t.table_name,\n" +
                "        NULL AS column_name,\n" +
                "        t.comment AS comment\n" +
                "    FROM system.information_schema.tables t\n" +
                "    WHERE t.table_schema = @schema\n" +
                "      AND t.table_type = 'BASE TABLE'\n" +
                ")\n" +
                "\n" +
                "UNION ALL\n" +
                "\n" +
                "-- VIEW comments\n" +
                "(\n" +
                "    SELECT \n" +
                "        v.table_catalog AS catalog_name,\n" +
                "        v.table_schema  AS table_schema,\n" +
                "        'view' AS object_type,\n" +
                "        v.table_name,\n" +
                "        NULL AS column_name,\n" +
                "        v.comment AS comment\n" +
                "    FROM system.information_schema.views v\n" +
                "    WHERE v.table_schema = @schema\n" +
                ")\n" +
                "\n" +
                "UNION ALL\n" +
                "\n" +
                "-- COLUMN comments (tables + views)\n" +
                "(\n" +
                "    SELECT \n" +
                "        c.table_catalog AS catalog_name,\n" +
                "        c.table_schema  AS table_schema,\n" +
                "        'column' AS object_type,\n" +
                "        c.table_name,\n" +
                "        c.column_name,\n" +
                "        c.comment AS comment\n" +
                "    FROM system.information_schema.columns c\n" +
                "    WHERE c.table_schema = @schema\n" +
                ")\n" +
                "\n" +
                "ORDER BY object_type, table_name, column_name;\n";
    }
}

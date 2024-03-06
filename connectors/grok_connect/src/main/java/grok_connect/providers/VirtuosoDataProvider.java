package grok_connect.providers;

import java.util.ArrayList;
import java.util.HashMap;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.utils.Property;
import serialization.Types;

public class VirtuosoDataProvider extends JdbcDataProvider {
    public VirtuosoDataProvider() {
        driverClassName = "virtuoso.jdbc4.Driver";

        descriptor = new DataSource();
        descriptor.type = "Virtuoso";
        descriptor.description = "Query Virtuoso database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.canBrowseSchema = true;
        descriptor.typesMap = new HashMap<String, String>() {{
            put("long varchar", Types.OBJECT);
            put("long nvarchar", Types.OBJECT);
            put("varchar", Types.STRING);
            put("nvarchar", Types.STRING);
            put("date", Types.DATE_TIME);
            put("datetime", Types.DATE_TIME);
            put("timestamp", Types.DATE_TIME);
            put("time", Types.DATE_TIME);
            put("integer", Types.INT);
            put("bigint", Types.BIG_INT);
            put("decimal", Types.FLOAT);
            put("real", Types.FLOAT);
            put("double precision", Types.FLOAT);
            put("xmltype", Types.OBJECT);
            put("any", Types.OBJECT);
            put("#.*binary", Types.BLOB);
            put("iri_id", Types.OBJECT);
        }};
    }

    @Override
    public String getConnectionString(DataConnection conn) {
        String connString = super.getConnectionString(conn);
        connString = connString.endsWith("/") ? connString : connString + "/";
        connString += "UID=" + conn.credentials.getLogin() + "/PWD=" + conn.credentials.getPassword();
        return connString;
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:virtuoso://" + conn.getServer() + port +
                "/TIMEOUT=100" + (conn.ssl() ? "/SSL" : "") + "/";
    }

    @Override
    public String getSchemasSql(String db) {
        return "SELECT DISTINCT table_schema FROM information_schema.columns";
    }

    @Override
    public String getSchemaSql(String db, String schema, String table) {
        String whereClause = String.format(" WHERE%s%s%s",
                db == null || db.isEmpty() ? "" : String.format(" LOWER(c.table_catalog) = LOWER('%s')", db),
                schema == null || schema.isEmpty() ? "" : String.format("%s c.table_schema = '%s'", db == null || db.isEmpty() ? "" : " AND",schema),
                table == null || table.isEmpty() ? "" : String.format("%s c.table_name = '%s'", (db == null || db.isEmpty())
                        && (schema == null || schema.isEmpty()) ? "" : " AND", table));
        return String.format("SELECT c.table_schema as table_schema, c.table_name as table_name, c.column_name as column_name, "
                        + "c.data_type as data_type FROM information_schema.columns c%s"
                , (db == null || db.isEmpty()) && schema == null && table == null ? "" : whereClause);
    }

    @Override
    public String limitToSql(String query, Integer limit) {
        return query + "top " + limit.toString() + " ";
    }

    @Override
    protected String getRegexQuery(String columnName, String regexExpression) {
        return String.format("(RDF_REGEX(%s, '%s') = 1)", columnName, regexExpression);
    }
}

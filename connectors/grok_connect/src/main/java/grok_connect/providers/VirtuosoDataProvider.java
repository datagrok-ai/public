package grok_connect.providers;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.Property;
import grok_connect.utils.ProviderManager;
import serialization.Types;

public class VirtuosoDataProvider extends JdbcDataProvider {
    public VirtuosoDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "virtuoso.jdbc4.Driver";

        descriptor = new DataSource();
        descriptor.type = "Virtuoso";
        descriptor.description = "Query Virtuoso database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.canBrowseSchema = true;
        descriptor.typesMap = new HashMap<String, String>() {{
            put("#.*char.*", Types.STRING);
            put("date", Types.DATE_TIME);
            put("#timestamp.*", Types.DATE_TIME);
            put("#time.*", Types.DATE_TIME);
            put("#interval.*", Types.STRING);
            put("int", Types.BIG_INT);
            put("#numeric.*", Types.FLOAT);
            put("float", Types.FLOAT);
            put("#geometry.*", Types.OBJECT);
            put("#geography.*", Types.OBJECT);
            put("uuid", Types.STRING);
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
        return query.replaceFirst("select", String.format("select top %s", limit));
    }

    @Override
    protected boolean isBigInt(int type, String typeName, int precision, int scale) {
        return type == java.sql.Types.OTHER && precision == 19 && scale == 0;
    }

    @Override
    protected String getRegexQuery(String columnName, String regexExpression) {
        return String.format("(RDF_REGEX(%s, '%s') = 1)", columnName, regexExpression);
    }
}

package grok_connect.providers;

import java.io.IOException;
import java.sql.SQLException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Properties;
import java.util.stream.Collectors;

import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.Property;
import grok_connect.utils.ProviderManager;
import grok_connect.utils.QueryCancelledByUser;
import serialization.Column;
import serialization.DataFrame;
import serialization.IntColumn;

public class VerticaDataProvider extends JdbcDataProvider {
    public VerticaDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "com.vertica.jdbc.Driver";

        descriptor = new DataSource();
        descriptor.type = "Vertica";
        descriptor.description = "Query Vertica database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.canBrowseSchema = true;
    }

    @Override
    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString() && conn.ssl())
            properties.setProperty("SSL", "true");
        return properties;
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:vertica://" + conn.getServer() + port + "/" + conn.getDb();
    }

    @Override
    public String getSchemasSql(String db) {
        return "SELECT DISTINCT schema_name as table_schema FROM v_catalog.schemata;";
    }

    @Override
    public String getSchemaSql(String db, String schema, String table) {
        String whereClause = String.format(" WHERE%s%s",
                schema == null ? "" : String.format(" c.table_schema = '%s'", schema),
                table == null ? "" : String.format("%s c.table_name = '%s'", schema == null ? "" : " AND", table));
        return String.format("SELECT c.table_schema as table_schema, c.table_name as table_name, c.column_name as column_name, \n"
                + "c.data_type as data_type, case t.table_type when 'VIEW' then 1 else 0 end as is_view \n"
                + "FROM v_catalog.columns c \n"
                + "JOIN v_catalog.all_tables t ON t.table_name = c.table_name%s UNION ALL "
                + "SELECT c.table_schema as table_schema, c.table_name as table_name, c.column_name as column_name, "
                + "c.data_type as data_type, case t.table_type when 'VIEW' then 1 else 0 end as is_view "
                + "FROM v_catalog.view_columns c JOIN v_catalog.all_tables t ON t.table_name = c.table_name%s;", schema == null && table == null ? "" : whereClause,
                schema == null && table == null ? "" : whereClause);
    }

    @Override
    protected String getRegexQuery(String columnName, String regexExpression) {
        return String.format("REGEXP_LIKE(%s, '%s')", columnName, regexExpression);
    }
}

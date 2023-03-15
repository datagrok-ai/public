package grok_connect.providers;

import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Properties;
import java.util.stream.Collectors;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataQuery;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncParam;
import grok_connect.table_query.AggrFunctionInfo;
import grok_connect.table_query.Stats;
import grok_connect.utils.Property;
import grok_connect.utils.ProviderManager;
import serialization.Types;

public class MySqlDataProvider extends JdbcDataProvider {
    public MySqlDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "com.mysql.cj.jdbc.Driver";

        descriptor = new DataSource();
        descriptor.type = "MySQL";
        descriptor.description = "Query MySQL database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.canBrowseSchema = true;
        descriptor.nameBrackets = "`";

        descriptor.typesMap = new HashMap<String, String>() {{
            put("bool", Types.BOOL);
            put("boolean", Types.BOOL);
            put("#bit.*", Types.INT);
            put("#tinyint.*", Types.INT);
            put("#smallint.*", Types.INT);
            put("#mediumint.*", Types.INT);
            put("#int.*", Types.INT);
            put("#integer.*", Types.INT);
            put("#bigint.*", Types.BIG_INT);
            put("#dec.*", Types.FLOAT);
            put("#decimal.*", Types.FLOAT);
            put("#float.*", Types.FLOAT);
            put("#double.*", Types.FLOAT);
            put("#double precision.*", Types.FLOAT);
            put("#char.*", Types.STRING);
            put("#varchar.*", Types.STRING);
            put("#text.*", Types.STRING);
            put("#mediumtext.*", Types.STRING);
            put("#longtext.*", Types.STRING);
            put("#tinytext.*", Types.STRING);
            put("date", Types.DATE_TIME);
            put("#datetime.*", Types.DATE_TIME);
            put("#timestamp.*", Types.DATE_TIME);
            put("#time.*", Types.DATE_TIME);
            put("year", Types.DATE_TIME);
        }};
        descriptor.aggregations.add(new AggrFunctionInfo(Stats.STDEV, "std(#)", Types.dataFrameNumericTypes));
    }

    @Override
    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString()) {
            properties.setProperty("zeroDateTimeBehavior", "convertToNull");
            if (conn.ssl()) {
                properties.setProperty("useSSL", "true");
                properties.setProperty("verifyServerCertificate", "false");
            }
        }
        return properties;
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:mysql://" + conn.getServer() + port + "/" + conn.getDb();
    }

    @Override
    public String getSchemasSql(String db) {
        return "SELECT DISTINCT table_schema FROM information_schema.columns";
    }

    @Override
    public String getSchemaSql(String db, String schema, String table) {
        List<String> filters = new ArrayList<>();

        if (db != null && db.length() != 0)
            filters.add("c.table_schema = '" + db + "'");

        if (table != null)
            filters.add("(c.table_name = '" + table + "')");

        String whereClause = filters.size() != 0 ? "WHERE " + String.join(" AND \n", filters) : "";

        return "SELECT c.table_schema as table_schema, c.table_name as table_name, c.column_name as column_name, "
                + "c.data_type as data_type, "
                + "case t.table_type when 'VIEW' then 1 else 0 end as is_view FROM information_schema.columns c "
                + "JOIN information_schema.tables t ON t.table_name = c.table_name " + whereClause +
                " ORDER BY c.table_name";
    }

    @Override
    protected String getRegexQuery(String columnName, String regexExpression) {
        return String.format("%s REGEXP '%s'", columnName, regexExpression);
    }

    @Override
    protected void appendQueryParam(DataQuery dataQuery, String paramName, StringBuilder queryBuffer) {
        FuncParam param = dataQuery.getParam(paramName);
        if (param.propertyType.equals("list")) {
            @SuppressWarnings("unchecked")
            List<String> values = ((ArrayList<String>) param.value);
            queryBuffer.append(values.stream().map(value -> "?").collect(Collectors.joining(", ")));
        } else {
            queryBuffer.append("?");
        }
    }

    @Override
    protected int setArrayParamValue(PreparedStatement statement, int n, FuncParam param) throws SQLException {
        @SuppressWarnings (value="unchecked")
        ArrayList<Object> lst = (ArrayList<Object>)param.value;
        if (lst == null || lst.size() == 0) {
            statement.setObject(n, null);
            return 0;
        }
        for (int i = 0; i < lst.size(); i++) {
            statement.setObject(n + i, lst.get(i));
        }
        return lst.size() - 1;
    }
}

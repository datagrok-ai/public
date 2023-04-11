package grok_connect.providers;

import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Types;
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
import grok_connect.resultset.ResultSetManager;
import grok_connect.table_query.AggrFunctionInfo;
import grok_connect.table_query.Stats;
import grok_connect.utils.Property;
import grok_connect.utils.ProviderManager;

public class Db2DataProvider extends JdbcDataProvider {
    public Db2DataProvider(ResultSetManager resultSetManager, ProviderManager providerManager) {
        super(resultSetManager, providerManager);
        driverClassName = "com.ibm.db2.jcc.DB2Driver";

        descriptor = new DataSource();
        descriptor.type = "DB2";
        descriptor.description = "Query DB2 database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.canBrowseSchema = true;
        descriptor.typesMap = new HashMap<String, String>() {{
            put("smallint", serialization.Types.INT);
            put("int", serialization.Types.INT);
            put("integer", serialization.Types.INT);
            put("bigint", serialization.Types.BIG_INT);
            put("real", serialization.Types.FLOAT);
            put("double", serialization.Types.FLOAT);
            put("decimal", serialization.Types.FLOAT);
            put("decfloat", serialization.Types.FLOAT);
            put("character", serialization.Types.STRING);
            put("varchar", serialization.Types.STRING);
            put("graphic", serialization.Types.STRING);
            put("vargraphic", serialization.Types.STRING);
            put("boolean", serialization.Types.BOOL);
            put("date", serialization.Types.DATE_TIME);
            put("timestamp", serialization.Types.DATE_TIME);
            put("time", serialization.Types.DATE_TIME);
            put("xml", serialization.Types.OBJECT);
            put("clob", serialization.Types.OBJECT);
            put("dbclob", serialization.Types.OBJECT);
        }};
        descriptor.aggregations.add(new AggrFunctionInfo(Stats.STDEV, "stddev(#)", serialization.Types.dataFrameNumericTypes));
    }

    @Override
    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString() && conn.ssl())
            properties.setProperty("sslConnection", "true");
        return properties;
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:db2://" + conn.getServer() + port + "/" + conn.getDb();
    }

    @Override
    public String getSchemasSql(String db) {
        return "SELECT DISTINCT TRIM(TABSCHEMA) AS table_schema FROM SYSCAT.COLUMNS;";
    }

    @Override
    public String getSchemaSql(String db, String schema, String table) {
        return String.format("SELECT C.TABSCHEMA as table_schema, C.TABNAME as table_name, C.COLNAME as column_name, "
                + "C.TYPENAME as data_type, case when T.TYPE = 'V' then 1 else 0 end as is_view FROM SYSCAT.COLUMNS C "
                + "JOIN SYSCAT.TABLES T ON C.TABNAME = T.TABNAME WHERE C.TABSCHEMA = '%s'%s ORDER BY C.COLNAME;",
                schema, table == null ? "" : String.format(" AND C.TABNAME = '%s'", table));
    }

    @Override
    protected String getRegexQuery(String columnName, String regexExpression) {
        return String.format("REGEXP_LIKE(%s, '%s')", columnName, regexExpression);
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

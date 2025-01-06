package grok_connect.providers;

import grok_connect.managers.ColumnManager;
import grok_connect.managers.bigint_column.ClickHouseBigIntColumnManager;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncParam;
import grok_connect.resultset.DefaultResultSetManager;
import grok_connect.resultset.ResultSetManager;
import grok_connect.utils.PatternMatcher;
import grok_connect.utils.PatternMatcherResult;
import grok_connect.utils.Property;
import serialization.Types;

import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Timestamp;
import java.util.*;

public class ClickHouseProvider extends JdbcDataProvider {
    public ClickHouseProvider() {
        driverClassName = "com.clickhouse.jdbc.ClickHouseDriver";
        descriptor = new DataSource();
        descriptor.type = "ClickHouse";
        descriptor.description = "Query ClickHouse database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.canBrowseSchema = true;
        descriptor.nameBrackets = "\"";
        descriptor.typesMap = new HashMap<String, String>() {{
            put("#^(int(8|16|32))$", serialization.Types.INT);
            put("#^(uint(8|16))$", serialization.Types.INT);
            put("#^(int(64|128|256))$", serialization.Types.BIG_INT);
            put("#^(uint(32|64|128|256))$", serialization.Types.BIG_INT);
            put("#float.*", serialization.Types.FLOAT);
            put("#decimal.*", serialization.Types.FLOAT);
            put("string", serialization.Types.STRING);
            put("uuid", serialization.Types.STRING);
            put("bool", serialization.Types.BOOL);
            put("#date.*", serialization.Types.DATE_TIME);
            put("#datetime.*", serialization.Types.DATE_TIME);
            put("#tuple.*", serialization.Types.OBJECT);
            put("#array.*", serialization.Types.OBJECT);
            put("#map.*", serialization.Types.OBJECT);
            put("point", serialization.Types.OBJECT);
            put("ring", serialization.Types.OBJECT);
            put("polygon", serialization.Types.OBJECT);
            put("multipolygon", serialization.Types.OBJECT);
        }};
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:clickhouse://" + conn.getServer() + port + "/" + conn.getDb();
    }


    @Override
    public String getSchemasSql(String db) {
        return "SELECT DISTINCT table_schema FROM information_schema.columns";
    }

    @Override
    public String getSchemaSql(String db, String schema, String table) {
        String whereClause = String.format(" WHERE%s%s",
                schema == null ? "" : String.format(" c.table_schema = '%s'", schema),
                table == null ? "" : String.format("%s c.table_name = '%s'", schema == null ? "" : " AND", table));
        return String.format("SELECT c.table_schema as table_schema, c.table_name as table_name, c.column_name as column_name, "
                        + "c.data_type as data_type, "
                        + "if(t.table_type ='VIEW', toInt8(1), toInt8(0)) as is_view FROM information_schema.columns c "
                        + "JOIN information_schema.tables t ON t.table_name = c.table_name%s"
                , schema == null && table == null ? "" : whereClause);
    }

    @Override
    protected String getRegexQuery(String columnName, String regexExpression) {
        return String.format("match(%s, '%s')", columnName, regexExpression);
    }

    @Override
    public PatternMatcherResult dateTimePatternConverter(String paramName, PatternMatcher matcher) {
        PatternMatcherResult result = new PatternMatcherResult();
        String queryFormat = "(toDateTime(%s) %s toDateTime(@%s))";
        String columnName = matcher.colName;
        switch (matcher.op) {
            case PatternMatcher.EQUALS:
                result.setQuery(String.format(queryFormat, columnName, "=", paramName));
                result.addParam(new FuncParam("datetime", paramName, matcher.values.get(0)));
                break;
            case PatternMatcher.BEFORE:
            case PatternMatcher.AFTER:
                result.setQuery(String.format(queryFormat, columnName,
                        PatternMatcher.cmp(matcher.op, matcher.include1), paramName));
                result.addParam(new FuncParam("datetime", paramName, matcher.values.get(0)));
                break;
            case PatternMatcher.RANGE_DATE_TIME:
                String name0 = paramName + "R0";
                String name1 = paramName + "R1";
                result.setQuery(String.format("(toDateTime(%s) %s toDateTime(@%s) "
                                + "AND toDateTime(%s) %s toDateTime(@%s))", columnName,
                        PatternMatcher.cmp(PatternMatcher.AFTER, matcher.include1),
                        name0, columnName, PatternMatcher.cmp(PatternMatcher.BEFORE, matcher.include2), name1));
                result.addParam(new FuncParam("datetime", name0, matcher.values.get(0)));
                result.addParam(new FuncParam("datetime", name1, matcher.values.get(1)));
                break;
            case PatternMatcher.NONE:
            default:
                result.setQuery("(1 = 1)");
                break;
        }
        return result;
    }

    @Override
    public ResultSetManager getResultSetManager() {
        Map<String, ColumnManager<?>> defaultManagersMap = DefaultResultSetManager.getDefaultManagersMap();
        defaultManagersMap.put(Types.BIG_INT, new ClickHouseBigIntColumnManager());
        return DefaultResultSetManager.fromManagersMap(defaultManagersMap);
    }

    @Override
    public void setDateTimeValue(FuncParam funcParam, PreparedStatement statement, int parameterIndex) throws SQLException {
        Calendar calendar = javax.xml.bind.DatatypeConverter.parseDateTime((String)funcParam.value);
        calendar.setTimeZone(TimeZone.getTimeZone("UTC"));
        Timestamp ts = new Timestamp(calendar.getTime().getTime());
        statement.setTimestamp(parameterIndex, ts);
    }
}

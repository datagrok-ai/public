package grok_connect.providers;

import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncParam;
import grok_connect.resultset.ResultSetManager;
import grok_connect.utils.PatternMatcher;
import grok_connect.utils.PatternMatcherResult;
import grok_connect.utils.Property;
import grok_connect.utils.ProviderManager;
import java.lang.reflect.Array;
import java.sql.Types;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class ClickHouseProvider extends JdbcDataProvider {
    public ClickHouseProvider(ResultSetManager resultSetManager, ProviderManager providerManager) {
        super(resultSetManager, providerManager);
        driverClassName = "com.clickhouse.jdbc.ClickHouseDriver";

        descriptor = new DataSource();
        descriptor.type = "ClickHouse";
        descriptor.description = "Query ClickHouse database";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.canBrowseSchema = true;
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
    public String getSchemaSql(String db, String schema, String table)
    {
        List<String> filters = new ArrayList<String>() {{
            add("c.table_schema = '" + ((schema != null) ? schema : descriptor.defaultSchema) + "'");
        }};

        if (db != null && db.length() != 0)
            filters.add("c.table_catalog = '" + db + "'");

        if (table != null)
            filters.add("c.table_name = '" + table + "'");

        String whereClause = "WHERE " + String.join(" AND \n", filters);

        return "SELECT c.table_schema as table_schema, c.table_name as table_name, c.column_name as column_name, "
                + "c.data_type as data_type, "
                + "if(t.table_type ='VIEW', toInt8(1), toInt8(0)) as is_view FROM information_schema.columns c "
                + "JOIN information_schema.tables t ON t.table_name = c.table_name " + whereClause +
                " ORDER BY c.table_name";
    }

    @Override
    protected String getRegexQuery(String columnName, String regexExpression) {
        return String.format("match(%s, '%s')", columnName, regexExpression);
    }

    @Override
    public PatternMatcherResult dateTimePatternConverter(FuncParam param, PatternMatcher matcher) {
        PatternMatcherResult result = new PatternMatcherResult();
        String queryFormat = "(toDateTime(%s) %s toDateTime(@%s))";
        String columnName = matcher.colName;
        switch (matcher.op) {
            case PatternMatcher.EQUALS:
                result.setQuery(String.format(queryFormat, columnName, "=", param.name));
                result.addParam(new FuncParam("datetime", param.name, matcher.values.get(0)));
                break;
            case PatternMatcher.BEFORE:
            case PatternMatcher.AFTER:
                result.setQuery(String.format(queryFormat, columnName,
                        PatternMatcher.cmp(matcher.op, matcher.include1), param.name));
                result.addParam(new FuncParam("datetime", param.name, matcher.values.get(0)));
                break;
            case PatternMatcher.RANGE_DATE_TIME:
                String name0 = param.name + "R0";
                String name1 = param.name + "R1";
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
    protected boolean isBigInt(int type, String typeName, int precision, int scale) {
        return type == Types.BIGINT
                || typeName.equalsIgnoreCase("UInt32")
                || typeName.equalsIgnoreCase("UInt64")
                || typeName.equalsIgnoreCase("UInt128")
                || typeName.equalsIgnoreCase("UInt256")
                || typeName.equalsIgnoreCase("Int64")
                || typeName.equalsIgnoreCase("Int128")
                || typeName.equalsIgnoreCase("Int256");
    }

    @Override
    protected boolean isDecimal(int type, String typeName, int scale) {
        return type == Types.DECIMAL;
    }

    @Override
    protected boolean isInteger(int type, String typeName, int precision, int scale) {
        return type == Types.INTEGER
                || type == Types.SMALLINT
                || type == Types.TINYINT
                || typeName.equalsIgnoreCase("UInt16")
                || typeName.equalsIgnoreCase("UInt8");
    }

    @Override
    protected Object convertArrayType(Object value) {
        if (value == null) {
            return Arrays.toString(new Object[]{});
        }
        return getStringArrayRepresentation(value);
    }

    private String getStringArrayRepresentation(Object array) {
        int length = Array.getLength(array);
        StringBuilder builder = new StringBuilder("[");
        for (int i = 0; i < length; i++) {
            Object o = Array.get(array, i);
            if (o.getClass().isArray()) {
                builder.append(getStringArrayRepresentation(o));
            } else {
                builder.append(o);
            }
            if (i != length - 1) {
                builder.append(", ");
            }
        }
        builder.append("]");
        return builder.toString();
    }
}

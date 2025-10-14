package grok_connect.connectors_info;

import java.util.*;
import java.util.regex.*;
import serialization.*;
import grok_connect.utils.*;
import grok_connect.table_query.*;


public abstract class DataProvider
{
    public static String CONN_AVAILABLE = "ok";
    public static String QUERY_COUNT = "#queryCount";
    public static String QUERY_MEMORY_LIMIT_MB = "#queryMemoryLimitMB";
    public static String QUERY_TIMEOUT_SEC = "#queryTimeoutSec";

    public DataSource descriptor;

    public String outputCsv;

    public abstract boolean autoInterpolation();

    public abstract DataFrame execute(FuncCall queryRun) throws QueryCancelledByUser, GrokConnectException;

    public DataFrame getSchemas(DataConnection dataConnection)
            throws QueryCancelledByUser, GrokConnectException {
        throw new UnsupportedOperationException();
    }

    public DataFrame getSchema(DataConnection dataConnection, String schema, String table)
            throws QueryCancelledByUser, GrokConnectException {
        throw new UnsupportedOperationException();
    }

    public abstract PatternMatcherResult numericPatternConverter(String paramName, String typeName, PatternMatcher matcher);
    public abstract PatternMatcherResult stringPatternConverter(String paramName, PatternMatcher matcher);
    public abstract PatternMatcherResult dateTimePatternConverter(String paramName, PatternMatcher matcher);
    public abstract PatternMatcherResult boolPatternConverter(String paramName, PatternMatcher matcher);
    public abstract PatternMatcherResult bigIntPatternConverter(String paramName, PatternMatcher matcher);

    @SuppressWarnings("unchecked")
    private PatternMatcherResult patternToQueryParam(FuncCall queryRun, FuncParam param, String colName) {
        String type = param.options.get("pattern");
        Map<String, Object> patterns = (Map<String, Object>) queryRun.options.get("patterns");
        PatternMatcher matcher = new PatternMatcher((Map<String, Object>)patterns.get(param.name), colName);
        switch (type) {
            case "num":
            case "double":
            case "int":
                return numericPatternConverter(param.name, param.options.get("pattern"), matcher);
            case "string":
                return stringPatternConverter(param.name, matcher);
            case "datetime":
                return dateTimePatternConverter(param.name, matcher);
            case "bool":
                return boolPatternConverter(param.name, matcher);
            case "bigint":
                return bigIntPatternConverter(param.name, matcher);
            default:
                throw new UnsupportedOperationException();
        }
    }

    public String convertPatternParamsToQueryParams(FuncCall queryRun, String query) {
        Pattern pattern = Pattern.compile("@(\\w+)\\((\\S+?)\\)");
        Matcher matcher = pattern.matcher(query);
        StringBuilder queryBuffer = new StringBuilder();
        int idx = 0;
        List<FuncParam> patternParams = new ArrayList<>();
        List<FuncParam> queryParams = new ArrayList<>();

        while (matcher.find()) {
            queryBuffer.append(query, idx, matcher.start());
            FuncParam param = queryRun.func.getParam(matcher.group(1));
            PatternMatcherResult pmr = patternToQueryParam(queryRun, param, matcher.group(2));
            patternParams.add(param);
            queryParams.addAll(pmr.params);
            queryBuffer.append(pmr.query);
            idx = matcher.end();
        }

        queryRun.func.removeParams(patternParams);
        queryRun.func.addParams(queryParams);

        queryBuffer.append(query.substring(idx));
        query = queryBuffer.toString();
        return query;
    }

    public abstract void testConnection(DataConnection conn) throws GrokConnectException;

    public String queryTableSql(TableQuery query) {
        throw new UnsupportedOperationException();
    }
}

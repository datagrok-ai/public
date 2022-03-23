package grok_connect.connectors_info;

import java.sql.*;
import java.text.*;
import java.util.*;
import java.util.regex.*;
import java.io.IOException;
import serialization.*;
import grok_connect.utils.*;
import grok_connect.providers.*;
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

    public abstract DataFrame execute(FuncCall queryRun)
            throws ClassNotFoundException, SQLException, ParseException, IOException, QueryCancelledByUser, GrokConnectException;

    public DataFrame getSchemas(DataConnection dataConnection)
            throws ClassNotFoundException, SQLException, ParseException, IOException, QueryCancelledByUser, GrokConnectException {
        throw new UnsupportedOperationException();
    }

    public DataFrame getSchema(DataConnection dataConnection, String schema, String table)
            throws ClassNotFoundException, SQLException, ParseException, IOException, QueryCancelledByUser, GrokConnectException {
        throw new UnsupportedOperationException();
    }

    public abstract PatternMatcherResult numericPatternConverter(FuncParam param, PatternMatcher matcher);
    public abstract PatternMatcherResult stringPatternConverter(FuncParam param, PatternMatcher matcher);
    public abstract PatternMatcherResult dateTimePatternConverter(FuncParam param, PatternMatcher matcher);

    private PatternMatcherResult patternToQueryParam(FuncCall queryRun, FuncParam param, String colName) {
        String type = param.options.get("pattern");
        Map patterns = (Map)queryRun.options.get("patterns");
        PatternMatcher matcher = new PatternMatcher((Map)patterns.get(param.name), colName);
        if (type.equals("num") || type.equals("double") || type.equals("int")) {
            return numericPatternConverter(param, matcher);
        } else if (type.equals("string")) {
            return stringPatternConverter(param, matcher);
        } else if (type.equals("datetime")) {
            return dateTimePatternConverter(param, matcher);
        } else {
            throw new UnsupportedOperationException();
        }
    }

    public String convertPatternParamsToQueryParams(FuncCall queryRun, String query) {
        Pattern pattern = Pattern.compile("@(\\w+)\\((\\S+)\\)");
        Matcher matcher = pattern.matcher(query);
        StringBuilder queryBuffer = new StringBuilder();
        int idx = 0;
        List<FuncParam> patternParams = new ArrayList<>();
        List<FuncParam> queryParams = new ArrayList<>();

        while (matcher.find()) {
            queryBuffer.append(query.substring(idx, matcher.start()));
            FuncParam param = queryRun.func.getParam(matcher.group(1));
            PatternMatcherResult pmr = patternToQueryParam(queryRun, param, matcher.group(2));
            patternParams.add(param);
            queryParams.addAll(pmr.params);

            queryBuffer.append(pmr.query);
            idx = matcher.end();
        }

        queryRun.func.removeParams(patternParams);
        queryRun.func.addParams(queryParams);

        queryBuffer.append(query.substring(idx, query.length()));
        query = queryBuffer.toString();
        return query;
    }

    public abstract String testConnection(DataConnection conn) throws ClassNotFoundException, SQLException;


    public DataFrame queryTable(DataConnection conn, TableQuery query)
            throws ClassNotFoundException, SQLException, ParseException, IOException, QueryCancelledByUser, GrokConnectException {
        throw new UnsupportedOperationException();
    }

    public String queryTableSql(DataConnection conn, TableQuery query) {
        throw new UnsupportedOperationException();
    }
}

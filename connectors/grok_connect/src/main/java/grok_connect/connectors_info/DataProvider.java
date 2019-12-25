package grok_connect.connectors_info;

import java.sql.*;
import java.text.*;
import java.util.*;
import java.util.regex.*;
import serialization.*;
import grok_connect.utils.*;
import grok_connect.providers.*;
import grok_connect.table_query.*;


public abstract class DataProvider
{
    public static String CONN_AVAILABLE = "ok";
    public static String QUERY_COUNT = "#queryCount";
    public static String QUERY_MEMORY_LIMIT_MB = "#queryMemoryLimitMB";

    public DataSource descriptor;

    public abstract boolean isParametrized();

    public abstract DataFrame execute(FuncCall queryRun)
            throws ClassNotFoundException, SQLException, ParseException;

    public DataFrame getSchema(DataConnection dataConnection, String schema, String table)
            throws ClassNotFoundException, SQLException, ParseException {
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
        Pattern pattern = Pattern.compile("@(\\w+)\\((\\w+)\\)");
        Matcher matcher = pattern.matcher(query);
        StringBuilder queryBuffer = new StringBuilder();
        int idx = 0;
        while (matcher.find()) {
            queryBuffer.append(query.substring(idx, matcher.start()));
            FuncParam param = queryRun.func.getParam(matcher.group(1));
            PatternMatcherResult pmr = patternToQueryParam(queryRun, param, matcher.group(2));
            queryRun.func.removeParam(param);
            queryRun.func.addParams(pmr.params);
            queryBuffer.append(pmr.query);
            idx = matcher.end();
        }
        queryBuffer.append(query.substring(idx, query.length()));
        query = queryBuffer.toString();
        return query;
    }

    public abstract String testConnection(DataConnection conn) throws ClassNotFoundException, SQLException;

    public static DataProvider getByName(String name) {
        DataProvider provider = Providers.get(0);

        for (ListIterator<DataProvider> providers = Providers.listIterator(); providers.hasNext(); ) {
            DataProvider tmp = providers.next();
            if (tmp.descriptor.type.equals(name)) {
                provider = tmp;
                break;
            }
        }

        return provider;
    }

    public static List<DataProvider> Providers = new ArrayList<DataProvider>() {{
        add(new AccessDataProvider());
        add(new AthenaDataProvider());
        add(new BigQueryDataProvider());
        add(new CassandraDataProvider());
        add(new Db2DataProvider());
        add(new FirebirdDataProvider());
        add(new HBaseDataProvider());
        add(new HiveDataProvider());
        add(new Hive2DataProvider());
        add(new MariaDbDataProvider());
        add(new MongoDbDataProvider());
        add(new MsSqlDataProvider());
        add(new MySqlDataProvider());
        add(new OracleDataProvider());
        add(new PostgresDataProvider());
        add(new RedshiftDataProvider());
        add(new SQLiteDataProvider());
        add(new TeradataDataProvider());
        add(new VerticaDataProvider());
    }};

    public static List<String> getAllProvidersTypes() {
        List<String> types = new ArrayList<>();

        for (DataProvider provider : Providers)
            types.add(provider.descriptor.type);

        return types;
    }

    public DataFrame queryTable(DataConnection conn, TableQuery query)
            throws ClassNotFoundException, SQLException, ParseException {
        throw new UnsupportedOperationException();
    }

    public String queryTableSql(DataConnection conn, TableQuery query) {
        throw new UnsupportedOperationException();
    }
}

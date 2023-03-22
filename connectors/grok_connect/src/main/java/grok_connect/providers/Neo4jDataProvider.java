package grok_connect.providers;

import java.sql.Types;
import java.util.ArrayList;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncParam;
import grok_connect.utils.PatternMatcher;
import grok_connect.utils.PatternMatcherResult;
import grok_connect.utils.Property;
import grok_connect.utils.ProviderManager;

public class Neo4jDataProvider extends JdbcDataProvider {
    public Neo4jDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "org.neo4j.jdbc.Driver";

        descriptor = new DataSource();
        descriptor.type = "Neo4j";
        descriptor.description = "Query Neo4j database";
        descriptor.commentStart = "//";
        descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.aggregations = null;
    }

    @Override
    public String getConnectionString(DataConnection conn) {
        String connString = super.getConnectionString(conn);
        if (!conn.hasCustomConnectionString() && !conn.ssl())
            connString += "?nossl";
        return connString;
    }

    @Override
    public void prepareProvider() throws ClassNotFoundException {
        super.prepareProvider();
        Class.forName("org.neo4j.jdbc.bolt.BoltDriver");
        Class.forName("org.neo4j.jdbc.boltrouting.BoltRoutingNeo4jDriver");
        Class.forName("org.neo4j.jdbc.http.HttpDriver");
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        String port = (conn.getPort() == null) ? "" : ":" + conn.getPort();
        return "jdbc:neo4j:bolt://" + conn.getServer() + port + "/" + conn.getDb();
    }

    @Override
    protected String getInQuery(PatternMatcher matcher, String names) {
        if (matcher.op.equals(PatternMatcher.IN)) {
            return String.format("(%s %s [%s])", matcher.colName, matcher.op, names);
        }
        return String.format("(NOT %s IN [%s])", matcher.colName, names);
    }

    @Override
    protected boolean isInteger(int type, String typeName, int precision, int scale) {
        // Neo4j only has one integer type, it's 64bit, so map to BigIntColumn
        return false;
    }

    @Override
    protected boolean isBigInt(int type, String typeName, int precision, int scale) {
        return typeName.equals("INTEGER") || type == Types.INTEGER;
    }

    @Override
    protected String getRegexQuery(String columnName, String regexExpression) {
        return String.format("%s =~ '%s'", columnName, regexExpression);
    }

    @Override
    public PatternMatcherResult stringPatternConverter(FuncParam param, PatternMatcher matcher) {
        PatternMatcherResult result = new PatternMatcherResult();
        String type = "string";
        String formatQuery = "(toLower(%s) %s toLower(@%s))";
        String value = ((String)matcher.values.get(0)).toLowerCase();
        if (matcher.op.equals(PatternMatcher.EQUALS)) {
            result.query = String.format(formatQuery, matcher.colName, "=", param.name);
            result.params.add(new FuncParam(type, param.name, value));
        } else if (matcher.op.equals(PatternMatcher.CONTAINS)) {
            result.query = String.format(formatQuery, matcher.colName, "CONTAINS", param.name);
            result.params.add(new FuncParam(type, param.name, value));
        } else if (matcher.op.equals(PatternMatcher.STARTS_WITH)) {
            result.query = String.format(formatQuery, matcher.colName, "STARTS WITH", param.name);
            result.params.add(new FuncParam(type, param.name, value));
        } else if (matcher.op.equals(PatternMatcher.ENDS_WITH)) {
            result.query = String.format(formatQuery, matcher.colName, "ENDS WITH", param.name);
            result.params.add(new FuncParam(type, param.name, value));
        } else if (matcher.op.equals(PatternMatcher.REGEXP)) {
            result.query = getRegexQuery(matcher.colName, value);
            result.params.add(new FuncParam(type, param.name, value));
        } else if (matcher.op.equals(PatternMatcher.IN) || matcher.op.equals(PatternMatcher.NOT_IN)) {
            String names = paramToNamesString(param, matcher, type, result);
            result.query = getInQuery(matcher, names);
        } else {
            result.query = "(1 = 1)";
        }
        return result;
    }
}

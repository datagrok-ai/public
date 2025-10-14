package grok_connect.table_query;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.DataQuery;
import grok_connect.connectors_info.FuncParam;
import grok_connect.providers.JdbcDataProvider;
import grok_connect.utils.GrokConnectUtil;
import grok_connect.utils.PatternMatcherResult;
import serialization.Types;

public class TableQuery extends DataQuery {
    public String whereOp = "and";
    public String havingOp = "and";
    public List<String> fields = new ArrayList<>();
    public List<String> pivots = new ArrayList<>();
    public List<FieldPredicate> whereClauses = new ArrayList<>();

    public List<GroupAggregation> aggregations = new ArrayList<>();
    public List<String> groupByFields = new ArrayList<>();
    public List<HavingPredicate> having = new ArrayList<>();
    public List<FieldOrder> orderBy = new ArrayList<>();
    List<TableJoin> joins = new ArrayList<>();
    public String tableName;
    public String schema;
    public Integer limit;

    public TableQuery() {
    }

    private String getFullTableName(String tableName, JdbcDataProvider provider) {
        tableName = provider.addBrackets(tableName);
        return GrokConnectUtil.isNotEmpty(schema) && !connection.dataSource.equals("SQLite")  && !connection.dataSource.equals("Databricks") && !tableName.contains(".")
                ? provider.addBrackets(schema) + "." + tableName : tableName;
    }

    public String toSql() {
        JdbcDataProvider provider = GrokConnect.providerManager.getByName(connection.dataSource);
        StringBuilder sql = new StringBuilder();
        StringBuilder sqlHeader = new StringBuilder();
        String table = tableName;
        if (table.contains(".")) {
            int idx = table.indexOf(".");
            schema = table.substring(0, idx);
            table = table.substring(idx + 1);
        }
        table = getFullTableName(table, provider);
        sql.append("SELECT");
        sql.append(System.lineSeparator());
        if (limit != null && !provider.descriptor.limitAtEnd) {
            sql.append(provider.limitToSql("", limit));
            sql.append(System.lineSeparator());
        }
        sql.append(getSelectFields(provider));
        sql.append("FROM");
        sql.append(System.lineSeparator());
        sql.append(table);

        if (!joins.isEmpty()) {
            sql.append(System.lineSeparator());
            for (TableJoin joinTable: joins) {
                sql.append(joinTable.joinType)
                        .append(" join ")
                        .append(getFullTableName(joinTable.rightTableName, provider));
                if (GrokConnectUtil.isNotEmpty(joinTable.rightTableAlias))
                    sql.append(" as ")
                            .append(provider.addBrackets(joinTable.rightTableAlias));
                sql.append(" on ");
                for (int i = 0; i < joinTable.leftTableKeys.size(); i++) {
                    if (i > 0) {
                        sql.append(" AND ");
                        sql.append(System.lineSeparator());
                    }
                    sql.append(getFullTableName(joinTable.leftTableName, provider))
                            .append('.')
                            .append(provider.addBrackets(joinTable.leftTableKeys.get(i)))
                            .append(" = ")
                            .append(GrokConnectUtil.isNotEmpty(joinTable.rightTableAlias) ? provider.addBrackets(joinTable.rightTableAlias) : getFullTableName(joinTable.rightTableName, provider))
                            .append(".")
                            .append(provider.addBrackets(joinTable.rightTableKeys.get(i)))
                            .append(System.lineSeparator());
                }
            }
        }

        if (!whereClauses.isEmpty()) {
            sql.append(System.lineSeparator());
            List<String> clauses = new ArrayList<>();
            sql.append("WHERE");
            sql.append(System.lineSeparator());
            for (FieldPredicate clause: whereClauses)
                clauses.add(String.format("  (%s)", preparePredicate(clause, sqlHeader, provider)));
            sql.append(String.join(String.format(" %s%s", whereOp, System.lineSeparator()), clauses));
        }

        if (!groupByFields.isEmpty()) {
            sql.append(System.lineSeparator());
            sql.append("GROUP BY");
            sql.append(System.lineSeparator());
            sql.append(
                    Stream.concat(
                            groupByFields.stream().map(provider::addBrackets),
                            !having.isEmpty() ? having.stream()
                                    .filter((predicate) -> GrokConnectUtil.isEmpty(predicate.aggType) && !groupByFields.contains(predicate.field))
                                    .map((predicate -> provider.addBrackets(predicate.field))) : Stream.empty()
                            )
                            .collect(Collectors.joining(", "))
            );
        }

        if (!pivots.isEmpty()) {
            sql.append(System.lineSeparator());
            sql.append("PIVOT ON");
            sql.append(System.lineSeparator());
            sql.append(String.join(", ", pivots));
        }

        if (!having.isEmpty()) {
            sql.append(System.lineSeparator());
            sql.append("HAVING");
            sql.append(System.lineSeparator());
            List<String> clauses = new ArrayList<>();
            for (FieldPredicate clause: having)
                clauses.add(String.format("\t(%s)", preparePredicate(clause, sqlHeader, provider)));
            sql.append(String.join(String.format(" %s%s", havingOp, System.lineSeparator()), clauses));
        }

        if (!orderBy.isEmpty()) {
            List<String> orders = new ArrayList<>();
            sql.append(System.lineSeparator());
            sql.append("ORDER BY");
            sql.append(System.lineSeparator());
            for (FieldOrder order: orderBy) {
                String orderField = provider.addBrackets(order.field);
                orders.add(String.format("%s%s", orderField, order.asc ? " asc" : " desc"));
            }
            sql.append(String.join(", ", orders));
        }
        String result;
        if (limit != null && provider.descriptor.limitAtEnd) {
            sql.append(System.lineSeparator());
            result = provider.limitToSql(sql.toString(), limit);
        }
        else
            result = sql.toString().trim();
        String header = sqlHeader.toString().trim();
        return header.isEmpty() ? result : header + System.lineSeparator() + result;
    }

    private String getSelectFields(JdbcDataProvider provider) {
        // use list to preserve order
        List<String> preparedFields = new ArrayList<>();
        Set<String> uniqueNames = new HashSet<>();
        for (String field : fields) {
            String[] splitField = field.split("\\.");
            String fieldName = splitField[splitField.length - 1];
            String bracket;
            if (uniqueNames.contains(fieldName) && splitField.length > 1 /* table alias in field */)
                bracket = provider.addBrackets(field) + " as " + provider.addBrackets(field.replaceAll("\\.", "_"));
            else
                bracket = provider.addBrackets(field);
            uniqueNames.add(fieldName);
            // fallback for old code
            int num = 1;
            while (true) {
                if (!preparedFields.contains(bracket)) {
                    preparedFields.add(bracket);
                    break;
                }
                bracket = String.format("%s AS %s", provider.addBrackets(field),
                        provider.addBrackets(String.format("%s(%d)", field, num++)));
            }
        }
        preparedFields.addAll(getAggFuncs().stream().map(provider::aggrToSql).filter(Objects::nonNull).collect(Collectors.toList()));
        return preparedFields.isEmpty() ? "*" + System.lineSeparator() : preparedFields.stream()
                .collect(Collectors.joining(String.format(",%s", System.lineSeparator()), "", System.lineSeparator()));
    }

    private List<GroupAggregation> getAggFuncs() {
        List<GroupAggregation> aggrs = new ArrayList<>();
        for (GroupAggregation aggregation: aggregations)
            if (!aggregation.aggType.equals(Stats.KEY))
                aggrs.add(aggregation);
        return aggrs;
    }

    private String preparePredicate(FieldPredicate clause, StringBuilder sqlHeader, JdbcDataProvider provider) {
        String paramName = clause.getParamName();
        clause.matcher.colName = provider.addBrackets(clause.field);
        if (clause instanceof HavingPredicate && GrokConnectUtil.isNotEmpty(((HavingPredicate) clause).aggType)) {
            provider.descriptor.aggregations.stream()
                    .filter((a) -> a.functionName.equals(((HavingPredicate) clause).aggType))
                    .findFirst()
                    .ifPresent(info -> clause.matcher.colName = info.dbFunctionName
                            .replaceAll("#", clause.matcher.colName));
        }
        PatternMatcherResult result;
        switch (clause.dataType) {
            case Types.NUM:
            case Types.FLOAT:
            case Types.INT:
                result = provider.numericPatternConverter(paramName, clause.dataType, clause.matcher);
                break;
            case Types.STRING:
                result = provider.stringPatternConverter(paramName, clause.matcher);
                break;
            case Types.DATE_TIME:
                result = provider.dateTimePatternConverter(paramName, clause.matcher);
                break;
            case Types.BOOL:
                result = provider.boolPatternConverter(paramName, clause.matcher);
                break;
            case Types.BIG_INT:
                result = provider.bigIntPatternConverter(paramName, clause.matcher);
                break;
            default:
                throw new UnsupportedOperationException(clause.dataType + " is not supported");
        }
        params.removeIf((p) -> p.name.equals(paramName));
        params.addAll(result.params);
        for (FuncParam p: result.params) {
            sqlHeader.append(String.format("--input: %s %s", p.propertyType, p.name));
            sqlHeader.append(System.lineSeparator());
        }

        return result.query;
    }
}

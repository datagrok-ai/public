package grok_connect.table_query;

import java.util.*;
import java.util.stream.Collectors;
import grok_connect.connectors_info.DataQuery;
import grok_connect.connectors_info.FuncParam;
import grok_connect.utils.GrokConnectUtil;
import serialization.Types;

public class TableQuery extends DataQuery {
    public String whereOp = "and";
    public String havingOp = "and";
    public List<String> fields = new ArrayList<>();
    public List<String> pivots = new ArrayList<>();
    public List<FieldPredicate> whereClauses = new ArrayList<>();

    public List<GroupAggregation> aggregations = new ArrayList<>();
    public List<String> groupByFields = new ArrayList<>();
    public List<FieldPredicate> having = new ArrayList<>();
    public List<FieldOrder> orderBy = new ArrayList<>();
    List<TableJoin> joins = new ArrayList<>();
    public String tableName;
    public String schema;
    public Integer limit;

    public TableQuery() {
    }

    public String toSql(AggrToSql aggrToSql, PatternToSql patternToSql, LimitToSql limitToSql, AddBrackets addBrackets,
                        boolean limitAtEnd) {
        StringBuilder sql = new StringBuilder();
        StringBuilder sqlHeader = new StringBuilder();
        String table = tableName;
        if (table.contains(".")) {
            int idx = table.indexOf(".");
            schema = table.substring(0, idx);
            table = table.substring(idx + 1);
        }
        table = addBrackets.convert(table);
        table = schema != null && !schema.isEmpty() && !connection.dataSource.equals("SQLite")  && !connection.dataSource.equals("Databricks")
                ? addBrackets.convert(schema) + "." + table : table;
        sql.append("SELECT");
        sql.append(System.lineSeparator());
        if (limit != null && !limitAtEnd) {
            sql.append(limitToSql.convert("", limit));
            sql.append(System.lineSeparator());
        }
        sql.append(getSelectFields(aggrToSql, addBrackets));
        sql.append("FROM");
        sql.append(System.lineSeparator());
        sql.append(table);

        if (!joins.isEmpty()) {
            sql.append(System.lineSeparator());
            for (TableJoin joinTable: joins) {
                sql.append(joinTable.joinType)
                        .append(" join ")
                        .append(addBrackets.convert(joinTable.rightTableName));
                if (GrokConnectUtil.isNotEmpty(joinTable.rightTableAlias))
                    sql.append(" as ")
                            .append(addBrackets.convert(joinTable.rightTableAlias));
                sql.append(" on ");
                for (int i = 0; i < joinTable.leftTableKeys.size(); i++) {
                    if (i > 0) {
                        sql.append(" AND ");
                        sql.append(System.lineSeparator());
                    }
                    sql.append(addBrackets.convert(joinTable.leftTableName))
                            .append('.')
                            .append(addBrackets.convert(joinTable.leftTableKeys.get(i)))
                            .append(" = ")
                            .append(addBrackets.convert(GrokConnectUtil.isNotEmpty(joinTable.rightTableAlias) ? joinTable.rightTableAlias : joinTable.rightTableName))
                            .append(".")
                            .append(addBrackets.convert(joinTable.rightTableKeys.get(i)))
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
                clauses.add(String.format("  (%s)", preparePredicate(clause, patternToSql, sqlHeader)));
            sql.append(String.join(String.format(" %s%s", whereOp, System.lineSeparator()), clauses));
        }

        if (!groupByFields.isEmpty()) {
            sql.append(System.lineSeparator());
            sql.append("GROUP BY");
            sql.append(System.lineSeparator());
            sql.append(
                    groupByFields.stream().map(addBrackets::convert).collect(Collectors.joining(", ")));
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
                clauses.add(String.format("\t(%s)",  preparePredicate(clause, patternToSql, sqlHeader)));
            sql.append(String.join(String.format(" %s%s", havingOp, System.lineSeparator()), clauses));
        }

        if (!orderBy.isEmpty()) {
            sql.append(System.lineSeparator());
            List<String> orders = new ArrayList<>();
            sql.append("ORDER BY");
            sql.append(System.lineSeparator());
            for (FieldOrder order: orderBy) {
                String orderField = connection.dataSource.equals("Access") ?
                        String.format("[%s]", order.field) : String.format("\"%s\"", order.field);
                orders.add(String.format("%s%s", orderField, order.asc ? " asc" : " desc"));
            }
            sql.append(String.join(", ", pad(orders)));
        }
        String result;
        if (limit != null && limitAtEnd) {
            sql.append(System.lineSeparator());
            result = limitToSql.convert(sql.toString(), limit);
        }
        else
            result = sql.toString().trim();
        String header = sqlHeader.toString().trim();
        return header.isEmpty() ? result : header + System.lineSeparator() + result;
    }

    private String getSelectFields(AggrToSql aggrToSql, AddBrackets addBrackets) {
        // use list to preserve order
        List<String> preparedFields = new ArrayList<>();
        Set<String> uniqueNames = new HashSet<>();
        for (String field : fields) {
            String[] splitField = field.split("\\.");
            String fieldName = splitField[splitField.length - 1];
            String bracket;
            if (uniqueNames.contains(fieldName) && splitField.length > 1 /* table alias in field */)
                bracket = addBrackets.convert(field) + " as " + addBrackets.convert(field.replaceAll("\\.", "_"));
            else
                bracket = addBrackets.convert(field);
            uniqueNames.add(fieldName);
            // fallback for old code
            int num = 1;
            while (true) {
                if (!preparedFields.contains(bracket)) {
                    preparedFields.add(bracket);
                    break;
                }
                bracket = String.format("%s AS %s", addBrackets.convert(field),
                        addBrackets.convert(String.format("%s(%d)", field, num++)));
            }
        }
        preparedFields.addAll(getAggFuncs().stream().map(aggrToSql::convert).filter(Objects::nonNull).collect(Collectors.toList()));
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

    private List<String> pad(List<String> strings) {
        strings.replaceAll(s -> "  " + s);
        return strings;
    }

    private String preparePredicate(FieldPredicate clause, PatternToSql patternToSql, StringBuilder sqlHeader) {
        String part = patternToSql.convert(clause);
        String dataType = part.contains("=") ? clause.dataType : Types.STRING;
        String paramName = clause.getParamName();
        sqlHeader.append(String.format("--input: %s %s", dataType, paramName));
        sqlHeader.append(System.lineSeparator());
        FuncParam param = new FuncParam();
        param.name = paramName;
        param.propertyType = dataType;
        param.options = new HashMap<>();
        param.options.put("default", clause.pattern);
        param.options.put("pattern",  part.contains("=") ? null : clause.dataType);
        params.add(param);
        return part;
    }
}

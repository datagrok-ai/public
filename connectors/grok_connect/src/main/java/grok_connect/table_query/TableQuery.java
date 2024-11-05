package grok_connect.table_query;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;
import grok_connect.connectors_info.DataConnection;

public class TableQuery {
    public String whereOp = "and";
    public String havingOp = "and";
    public List<String> fields = new ArrayList<>();
    public List<String> pivots = new ArrayList<>();
    public List<FieldPredicate> whereClauses = new ArrayList<>();

    public List<GroupAggregation> aggregations = new ArrayList<>();
    public List<String> groupByFields = new ArrayList<>();
    public List<FieldPredicate> having = new ArrayList<>();
    public List<FieldOrder> orderBy = new ArrayList<>();
    public String tableName;
    public String schema;
    public DataConnection connection;
    public Integer limit;

    public TableQuery() {
    }

    public String toSql(AggrToSql aggrToSql, PatternToSql patternToSql, LimitToSql limitToSql, AddBrackets addBrackets,
                        boolean limitAtEnd) {
        StringBuilder sql = new StringBuilder();
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
        if (!whereClauses.isEmpty()) {
            sql.append(System.lineSeparator());
            List<String> clauses = new ArrayList<>();
            sql.append("WHERE");
            sql.append(System.lineSeparator());
            for (FieldPredicate clause: whereClauses)
                clauses.add(String.format("  (%s)", patternToSql.convert(clause)));
            sql.append(clauses.stream()
                    .collect(Collectors
                            .joining(String.format(" %s%s", whereOp, System.lineSeparator())))
            );
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
                clauses.add(String.format("\t(%s)", patternToSql.convert(clause)));
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
            result = sql.toString();
        return result;
    }

    private String getSelectFields(AggrToSql aggrToSql, AddBrackets addBrackets) {
        List<String> preparedFields = new ArrayList<>();
        for (String field : fields) {
            int num = 1;
            String bracket = addBrackets.convert(field);
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
        return preparedFields.isEmpty() ? "*\n" : preparedFields.stream()
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
}

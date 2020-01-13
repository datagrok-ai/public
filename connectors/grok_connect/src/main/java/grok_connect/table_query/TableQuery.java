package grok_connect.table_query;

import java.util.*;
import org.apache.commons.lang.StringUtils;


public class TableQuery
{
    public String tableName;
    public String schema;
    public List<String> fields = new ArrayList<>();
    public List<String> pivots = new ArrayList<>();

    /// Logical operator applied to all [where] conditions.
    /// Allowed values are 'or' and 'and'.
    public String whereOp = "and";

    /// Conditions applied to the rows before they are aggregated.
    /// Translates into SQL's WHERE clause.
    public List<FieldPredicate> whereClauses = new ArrayList<>();

    /// Aggregations to apply.
    public List<GroupAggregation> aggregations = new ArrayList<>();

    public List<String> groupByFields = new ArrayList<>();

    public List<GroupAggregation> getAggFuncs() {
        List<GroupAggregation> aggrs = new ArrayList<>();
        for (GroupAggregation aggregation: aggregations)
            if (!aggregation.aggType.equals(Stats.KEY))
                aggrs.add(aggregation);
        return aggrs;
    }

    /// Logical operator applied to all [having] conditions.
    /// Allowed values are 'or' and 'and'.
    public String havingOp = "and";

    /// Conditions applied to the rows before they are aggregated.
    /// Translates into SQL's HAVING clause.
    public List<FieldPredicate> having = new ArrayList<>();

    /// Order to apply.
    public List<FieldOrder> orderBy = new ArrayList<>();

    /// Number of records to return, null means there is no limit.
    public Integer limit;

    public TableQuery() {
    }

    public String toSql(AggrToSql aggrToSql, PatternToSql patternToSql, LimitToSql limitToSql, String nameBrackets,
                        boolean limitAtEnd) {
        if (aggrToSql == null)
            aggrToSql = GroupAggregation::toSqlString;
        if (patternToSql == null)
            patternToSql = (fp) -> fp.field + " " + fp.pattern;
        if (limitToSql == null)
            limitToSql = (query, limit) -> query + "limit " + limit.toString();

        String sql = "";
        String str = "";
        List<String> selectFields = new ArrayList<>(fields);
        for (GroupAggregation func : getAggFuncs())
            selectFields.add(aggrToSql.convert(func));
        if (selectFields.size() == 0)
            return "";
        selectFields = pad(selectFields);
        str += StringUtils.join(selectFields,",\n") + "\n";
        String _tableName = tableName;
        if (_tableName.contains(".")) {
            int idx = _tableName.indexOf(".");
            schema = _tableName.substring(0, idx);
            _tableName = _tableName.substring(idx + 1);
        }
        if (_tableName.contains(" "))
            _tableName = nameBrackets.substring(0, 1) + _tableName +
                    nameBrackets.substring(nameBrackets.length() - 1, nameBrackets.length());
        _tableName = (schema != null && schema.length() != 0) ? schema + "." + _tableName : _tableName;
        sql += "select \n" + ((limit != null && !limitAtEnd) ? limitToSql.convert("", limit) + "\n" : "") +
                str + "from \n  " + _tableName + "\n";

        if (!whereClauses.isEmpty()) {
            List<String> clauses = new ArrayList<>();
            sql += "where\n";
            for (FieldPredicate clause: whereClauses)
                clauses.add("  (" + patternToSql.convert(clause) + ")");
            sql += StringUtils.join(clauses, " " + whereOp + "\n") + "\n";
        }

        if (!groupByFields.isEmpty())
            sql += "group by\n  " + StringUtils.join(groupByFields, ", ") + "\n";

        // the actual pivoting is done on the client, this code is only for the query to look good
        if (!pivots.isEmpty())
            sql += "pivot on\n  " + StringUtils.join(pivots, ", ") + "\n";

        if (!having.isEmpty()) {
            sql += "having\n";
            sql += "";
        }

        if (!orderBy.isEmpty()) {
            List<String> orders = new ArrayList<>();
            sql += "order by\n";
            for (FieldOrder order: orderBy)
                orders.add(order.field + (order.asc ? " asc" : " desc"));
            sql += StringUtils.join(pad(orders), ", ") + "\n";
        }

        if (limit != null && limitAtEnd)
            sql = limitToSql.convert(sql, limit);

        return sql;
    }

    private static List<String> pad(List<String> strings) {
        for (int n = 0; n < strings.size(); n++)
            strings.set(n, "  " + strings.get(n));
        return strings;
    }
}

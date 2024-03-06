package grok_connect.utils;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class PatternMatcher {
    public static final String RANGE_NUM = "-";
    public static final String NONE = "none";
    public static final String CONTAINS = "contains";
    public static final String STARTS_WITH = "starts with";
    public static final String ENDS_WITH = "ends with";
    public static final String EQUALS = "equals";
    public static final String REGEXP = "regex";
    public static final String IN = "in";
    public static final String NOT_IN = "not in";
    public static final String BEFORE = "before";
    public static final String AFTER = "after";
    public static final String RANGE_DATE_TIME = "range";
    public static final String IS_NULL = "is null";
    public static final String IS_NOT_NULL = "is not null";
    /// Expression as entered by user.
    public String expression;
    /// Operation (such as "EQUALS", "BEFORE", etc). [commonOption] contains list of available options.
    public String op;
    /// Column name this expression should be applied to; can be empty (for instance when used in search patterns).
    public String colName;
    public List<String> values;
    public Boolean include1;
    public Boolean include2;

    @SuppressWarnings({"unchecked", "raw"})
    public PatternMatcher(Map<String, Object> matcher, String colName) {
        this.colName = colName;
        expression = (String)matcher.get("expression");
        op = (String)matcher.get("op");
        values = (List<String>) getOptionalValue(matcher, "values");
        include1 = (Boolean)getOptionalValue(matcher, "include1");
        include2 = (Boolean)getOptionalValue(matcher, "include2");
    }

    private Object getOptionalValue(Map<String, Object> matcher, String key) {
        return matcher.getOrDefault(key, null);
    }

    public String toSql(String type, String variable) {
        if (op.equals(NONE))
            return "(1 = 1)";
        switch (type) {
            case "num":
            case "double":
            case "int":
                switch (op) {
                    case RANGE_NUM:
                        return "(" + variable + " >= " + values.get(0) + " AND " + variable + " <= " + values.get(1) + ")";
                    case IN:
                        return "(" + variable + " IN " + "(" + String.join(",", values) + "))";
                    case NOT_IN:
                        return "(" + variable + " NOT IN " + "(" + String.join(",", values) + "))";
                }
                return "(" + variable + op + values.get(0) + ")";
            case "string":
                String val = values.get(0).toLowerCase();
                switch (op) {
                    case EQUALS:
                        return "(LOWER(" + variable + ") LIKE '" + val + "')";
                    case CONTAINS:
                        return "(LOWER(" + variable + ") LIKE '%" + val + "%')";
                    case STARTS_WITH:
                        return "(LOWER(" + variable + ") LIKE '" + val + "%')";
                    case ENDS_WITH:
                        return "(LOWER(" + variable + ") LIKE '%" + val + "')";
                    case REGEXP:
                        return "(LOWER(" + variable + ") SIMILAR TO '" + val + "')";
                    case IN:
                        List<String> _values = new ArrayList<>();
                        for (String v : values)
                            _values.add("'" + v.toLowerCase() + "'");
                        return "(LOWER(" + variable + ") IN (" + String.join(",", _values) + "))";
                    default:
                        throw new UnsupportedOperationException("Unknown operation " + op);
                }
            case "datetime":
                switch (op) {
                    case EQUALS:
                        return "(" + variable + " = '" + values.get(0) + "')";
                    case BEFORE:
                        return "(" + variable + cmp(BEFORE, include1) + "'" + values.get(0) + "')";
                    case AFTER:
                        return "(" + variable + cmp(AFTER, include1) + "'" + values.get(0) + "')";
                    case RANGE_DATE_TIME:
                        return "((" + variable + cmp(AFTER, include1) + "'" + values.get(0) + "') AND (" +
                                variable + cmp(BEFORE, include1) + "'" + values.get(0) + "'))";
                    default:
                        throw new UnsupportedOperationException("Unknown operation " + op);
                }
            default:
                throw new UnsupportedOperationException();
        }
    }

    public static String cmp(String op, boolean include) {
        if (op.equals(PatternMatcher.BEFORE)) {
            return include ? " <= " : " < ";
        } else if (op.equals(PatternMatcher.AFTER)) {
            return include ? " >= " : " > ";
        } else {
            return "";
        }
    }
}

package grok_connect.table_query;

public interface PatternToSql {
    public String convert(FieldPredicate fp);
}

package grok_connect.table_mutation;

/** The db-schemas {path, code, message} refusal shape + a live pre-check count. */
public class DestructiveAction {
    public String path; // "table" or "columns.<name>"
    public String code; // drop-table|truncate-table|drop-column|type-change|not-null-violation|required-without-default
    public String message;
    public Long count; // live rows / non-null values affected
}

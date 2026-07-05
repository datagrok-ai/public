package grok_connect.table_mutation;

/** Per-row error, byte-compatible with the domain-schemas error shape. */
public class RowError {
    public int index;
    public String column;
    public String code;
    public String message;
}

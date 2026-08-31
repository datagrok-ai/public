package grok_connect.table_mutation;

/** Column definition for CreateTable / AlterTable addColumn; plain Gson class, no #type. */
public class ColumnSpec {
    public String name;
    public String type; // dg type: string|int|bigint|float|bool|datetime
    public boolean nullable = true;
    public String defaultValue; // always a String on the wire; typed by {@link #type} at emission
}

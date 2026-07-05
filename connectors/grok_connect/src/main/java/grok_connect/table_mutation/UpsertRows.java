package grok_connect.table_mutation;

/** Insert-or-update; {@link #matchKeys} is required non-empty. */
public class UpsertRows extends InsertRows {
    public UpsertRows() {
        type = "UpsertRows";
    }
}

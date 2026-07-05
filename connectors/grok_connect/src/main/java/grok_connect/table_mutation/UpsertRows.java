package grok_connect.table_mutation;

/** Insert-or-update; {@link #matchKeys} is required non-empty. */
public class UpsertRows extends InsertRows {
    public UpsertRows() {
        type = "UpsertRows";
        mode = "upsert";
    }

    /** Promotes a plain {@link InsertRows} carrying {@code mode == "upsert"} into an upsert view. */
    public UpsertRows(InsertRows src) {
        this();
        tableName = src.tableName;
        schema = src.schema;
        catalog = src.catalog;
        connection = src.connection;
        columns = src.columns;
        columnTypes = src.columnTypes;
        rows = src.rows;
        matchKeys = src.matchKeys;
        keyColumns = src.keyColumns;
        allOrNothing = src.allOrNothing;
        errorOnDuplicate = src.errorOnDuplicate;
        returnGeneratedKeys = src.returnGeneratedKeys;
    }
}

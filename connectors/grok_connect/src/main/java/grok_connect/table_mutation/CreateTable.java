package grok_connect.table_mutation;

import java.util.List;

/** CREATE TABLE with optional primary key and indexes; FKs are added afterwards via {@link AddKey}. */
public class CreateTable extends DdlMutation {
    public List<ColumnSpec> columns;
    public List<String> primaryKey;
    public List<IndexSpec> indexes;
    public boolean ifNotExists = false;

    public CreateTable() {
        type = "CreateTable";
    }
}

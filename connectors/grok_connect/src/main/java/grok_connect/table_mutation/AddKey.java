package grok_connect.table_mutation;

import java.util.List;

/** ALTER TABLE ADD CONSTRAINT — primary or foreign key. */
public class AddKey extends DdlMutation {
    public String keyType; // primary|foreign
    public String constraintName; // auto-generated (pk_<table> / fk_<table>_<col>) when absent
    public List<String> columns;
    public String refTable;         // foreign only
    public List<String> refColumns; // foreign only

    public AddKey() {
        type = "AddKey";
    }
}

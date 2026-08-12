package grok_connect.table_mutation;

/** Single-action ALTER TABLE; per-action field presence is validated at emission time. */
public class AlterTable extends DdlMutation {
    public String action; // addColumn|dropColumn|renameColumn|changeType|setNullable
    public String columnName;
    public String newName;    // renameColumn
    public ColumnSpec column; // addColumn
    public String newType;    // changeType, a dg type
    public Boolean nullable;  // setNullable — an explicit value is required

    public AlterTable() {
        type = "AlterTable";
    }
}

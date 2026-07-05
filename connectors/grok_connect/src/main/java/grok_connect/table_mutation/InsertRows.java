package grok_connect.table_mutation;

import java.util.List;

public class InsertRows extends TableMutation {
    public List<String> columns;
    public List<String> columnTypes; // dg types per column
    public List<List<Object>> rows;  // rows == null && bulk: payload streams separately (WO-5)
    public boolean bulk;
    public String mode = "insert"; // insert|upsert|update (WO-6)
    public List<String> matchKeys;
    public List<String> keyColumns;
    public boolean allOrNothing = true;
    public boolean errorOnDuplicate;
    public boolean returnGeneratedKeys;

    public InsertRows() {
        type = "InsertRows";
    }
}

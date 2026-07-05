package grok_connect.table_mutation;

import java.util.List;
import java.util.Map;

import grok_connect.table_query.FieldPredicate;

public class UpdateRows extends TableMutation {
    public Map<String, Object> set;
    public List<String> setColumns; // SET column order; setTypes is aligned to it, not to map iteration order
    public List<String> setTypes;   // dg types per setColumns entry
    public List<FieldPredicate> whereClauses;
    public String whereOp = "and";

    public UpdateRows() {
        type = "UpdateRows";
    }
}

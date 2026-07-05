package grok_connect.table_mutation;

import java.util.List;

import grok_connect.table_query.FieldPredicate;

public class UpdateRows extends TableMutation {
    // parallel lists aligned to setColumns — JSON map ordering is not a contract across Gson/Dart
    public List<String> setColumns;
    public List<Object> setValues;
    public List<String> setTypes; // dg types
    public List<FieldPredicate> whereClauses;
    public String whereOp = "and";

    public UpdateRows() {
        type = "UpdateRows";
    }
}

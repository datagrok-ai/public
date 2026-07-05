package grok_connect.table_mutation;

import java.util.List;

import grok_connect.table_query.FieldPredicate;

public class DeleteRows extends TableMutation {
    public List<FieldPredicate> whereClauses;
    public String whereOp = "and";
    public boolean allowFullTable = false; // empty WHERE without this flag is a structured error

    public DeleteRows() {
        type = "DeleteRows";
    }
}

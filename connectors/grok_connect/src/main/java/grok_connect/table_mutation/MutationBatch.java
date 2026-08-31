package grok_connect.table_mutation;

import java.util.List;

/** Ordered operations executed in a single transaction. */
public class MutationBatch extends TableMutation {
    public List<TableMutation> operations;

    public MutationBatch() {
        type = "MutationBatch";
    }
}

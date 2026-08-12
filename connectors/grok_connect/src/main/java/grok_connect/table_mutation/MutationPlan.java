package grok_connect.table_mutation;

import java.util.List;

/** What a mutation would execute: exact per-provider SQL + destructive-action pre-check results. */
public class MutationPlan {
    public List<String> statements;
    public List<DestructiveAction> destructive;
    public Boolean transactionalDdl; // null for DML plans; from the descriptor for DDL
}

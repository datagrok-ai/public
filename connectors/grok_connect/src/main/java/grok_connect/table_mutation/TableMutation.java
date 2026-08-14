package grok_connect.table_mutation;

import grok_connect.connectors_info.DataQuery;

/** Base class for structured write operations; mirrors table_query.TableQuery addressing. */
public abstract class TableMutation extends DataQuery {
    public String tableName;
    public String schema;
    public String catalog;
    public boolean dryRun = false; // return MutationResult.plan without executing
    public boolean confirmDestructive = false; // required when the recomputed plan has destructive actions
}

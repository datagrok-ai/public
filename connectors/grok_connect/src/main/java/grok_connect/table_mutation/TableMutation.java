package grok_connect.table_mutation;

import grok_connect.connectors_info.DataQuery;

/** Base class for structured write operations; mirrors table_query.TableQuery addressing. */
public abstract class TableMutation extends DataQuery {
    public String tableName;
    public String schema;
    public String catalog;
}

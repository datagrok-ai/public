package grok_connect.table_query;


public class GroupAggregation
{
    /// Aggregation type. See [Stats].
    public String aggType;

    /// Name of the source column to aggregate
    public String colName;

    /// Name of the result (aggregated) column.
    public String resultColName;

    /// auxiliary field for workspace tree-originated drag-n-drop.
    public String srcTableName;

    public String getDstColName() {
        return resultColName != null ? resultColName : (aggType + "(" + colName + ")");
    }

    public GroupAggregation(String aggType, String colName, String resultColName) {
        this.aggType = aggType;
        this.colName = colName;
        this.resultColName = resultColName;
    }

    public String toSqlString() {
        String func = aggType.equals(Stats.TOTAL_COUNT) ? "count(*)" : getDstColName();
        return func + ((resultColName == null) ? "" : "as " + resultColName);
    }

    public String toString() {
        return getDstColName();
    }

    public String getName() {
        return aggType + "(" + colName + ")";
    }
}

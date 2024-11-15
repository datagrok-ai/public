package grok_connect.table_query;

import java.util.List;

public class TableJoin {
    public String leftTableName;
    public String rightTableName;
    public String rightTableAlias;
    public String joinType;
    public List<String> leftTableKeys;
    public List<String> rightTableKeys;

    public TableJoin () {
    }
}

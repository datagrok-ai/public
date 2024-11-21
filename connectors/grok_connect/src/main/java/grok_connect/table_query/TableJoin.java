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

    public TableJoin(String leftTableName, String rightTableName, String rightTableAlias, String joinType, List<String> leftTableKeys, List<String> rightTableKeys) {
        this.leftTableName = leftTableName;
        this.rightTableName = rightTableName;
        this.rightTableAlias = rightTableAlias;
        this.joinType = joinType;
        this.leftTableKeys = leftTableKeys;
        this.rightTableKeys = rightTableKeys;
    }
}

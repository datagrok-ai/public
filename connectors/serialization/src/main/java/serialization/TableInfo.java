package serialization;

import java.util.*;
import com.google.gson.annotations.SerializedName;


/// Represents [DataFrame] metadata.
public class TableInfo
{
    @SerializedName("#type")
    public String type = "TableInfo";

    public String id;
    public int[] datasets;
    public String name;

    public ColumnInfo[] columns;

    public Map<String, String> tags = new HashMap<>();
    public String queryRun;

    public int rowCount;
    public int colCount;

    public TableInfo(DataFrame t)
    {
        name = t.name;
        if (t.columns.size() > 0) {
            rowCount = t.columns.get(0).length;
            colCount = t.columns.size();
            columns = new ColumnInfo[colCount];
            for (int n = 0; n < colCount; n++)
                columns[n] = new ColumnInfo(t.columns.get(n));
        } else {
            rowCount = 0;
            colCount = 0;
        }
        tags = t.tags;
    }
}

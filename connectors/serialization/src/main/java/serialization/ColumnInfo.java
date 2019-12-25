package serialization;

import java.util.*;
import com.google.gson.annotations.SerializedName;


/// Represents [Column] metadata.
public class ColumnInfo {
    @SerializedName("#type")
    public String _type;

    public String id;
    public TableInfo tableInfo;
    public String name;
    public String type;

    public ColumnInfo(Column c) {
        _type = "ColumnInfo";
        name = c.name;
        type = c.getType();
    }
}

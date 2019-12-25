package grok_connect.table_query;

import java.util.*;
import com.google.gson.annotations.*;


public class AggrFunctionInfo
{
    @SerializedName("#type")
    public String _type = "AggrFunctionInfo";

    /// Whether to show it automatically in the property panel
    public boolean auto = true;

    /// Function alias (see [Stats])
    public String functionName;

    /// Defines the way the function should be rendered for SQL, using '#' as a column name.
    /// Examples:
    ///   Maximum:  'max(#)'
    ///   Number of empty values: 'count(*) - count(#)'
    public String dbFunctionName;

    /// List of types that the function is applicable to (see [Types])
    public List<String> srcTypes;

    /// Result type (see [Types])
    public String resType;

    public AggrFunctionInfo(String functionName, String dbFunctionName, List<String> srcTypes, String resType) {
        this.functionName = functionName;
        this.dbFunctionName = dbFunctionName;
        this.srcTypes = srcTypes;
        this.resType = resType;
    }

    public AggrFunctionInfo(String functionName, String dbFunctionName, List<String> srcTypes) {
        this.functionName = functionName;
        this.dbFunctionName = dbFunctionName;
        this.srcTypes = srcTypes;
        this.resType = resType;
    }

    public String toString() {
        return functionName;
    }
}

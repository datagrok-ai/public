package grok_connect.connectors_info;

import java.util.*;
import com.google.gson.annotations.*;
import grok_connect.utils.*;
import grok_connect.table_query.*;


public class DataSource
{
    @SerializedName("#type")
    public String _type = "DataSource";

    public String type;
    public String category = "Database";
    public String description;
    public String commentStart = "--";
    public String queryLanguage = "sql";
    public boolean canBrowseSchema = false;
    public String nameBrackets = "[]";
    public boolean limitAtEnd = true;

    public List<Property> connectionTemplate;
    public List<Property> credentialsTemplate;
    public List<Property> queryTemplate;

    public List<AggrFunctionInfo> aggregations;
    public Map<String, String> typesMap;
}

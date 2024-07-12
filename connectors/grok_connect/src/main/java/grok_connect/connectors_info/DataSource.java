package grok_connect.connectors_info;

import java.util.*;
import com.google.gson.annotations.*;
import grok_connect.utils.*;
import grok_connect.table_query.*;
import serialization.Types;


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
    public String defaultSchema = null;
    public String nameBrackets = "[]";
    public boolean limitAtEnd = true;

    public List<Property> connectionTemplate;
    public List<Property> credentialsTemplate;
    public List<Property> cacheTemplate = new ArrayList<Property>() {{
        add(new Property(Property.BOOL_TYPE, DbCredentials.CACHE_SCHEMA, "Cache results of schema browsing"));
        add(new Property(Property.BOOL_TYPE, DbCredentials.CACHE_RESULTS, "Cache query results"));
        add(new Property(Property.STRING_TYPE, DbCredentials.CACHE_INVALIDATE_SCHEDULE,
                "Invalidation schedule for cache. Use cron expression to set valid invalidation schedule"));
    }};

    public List<Property> queryTemplate;

    public List<AggrFunctionInfo> aggregations = new ArrayList<AggrFunctionInfo>() {{
        add(new AggrFunctionInfo(Stats.AVG, "avg(#)", Types.dataFrameNumericTypes));
        add(new AggrFunctionInfo(Stats.MIN, "min(#)", Types.dataFrameNumericTypes));
        add(new AggrFunctionInfo(Stats.MAX, "max(#)", Types.dataFrameNumericTypes));
        add(new AggrFunctionInfo(Stats.SUM, "sum(#)", Types.dataFrameNumericTypes));
        add(new AggrFunctionInfo(Stats.TOTAL_COUNT, "count(*)", Types.dataFrameColumnTypes));
        add(new AggrFunctionInfo(Stats.VALUE_COUNT, "count(#)", Types.dataFrameColumnTypes));
        add(new AggrFunctionInfo(Stats.MISSING_VALUE_COUNT, "count(*) - count(#)", Types.dataFrameColumnTypes));
    }};

    public Map<String, String> typesMap;
}

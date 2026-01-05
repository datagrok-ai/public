package grok_connect.providers.utils;

import com.google.gson.internal.LinkedTreeMap;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataQuery;
import grok_connect.connectors_info.FuncCall;
import grok_connect.connectors_info.FuncParam;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Builder fore reducing amount of code
 */
public class FuncCallBuilder {
    private static FuncCallBuilder builder;
    private FuncCall funcCall;

    private FuncCallBuilder() {
    }

    public static FuncCallBuilder getBuilder() {
        if (builder == null) {
            builder = new FuncCallBuilder();
            builder.restore();
        }
        return builder;
    }

    public static FuncCall fromQuery(String query) {
        return getBuilder().addQuery(query).build();
    }

    public void restore() {
        funcCall = new FuncCall();
        funcCall.func = new DataQuery();
        funcCall.func.params = new ArrayList<>();
        funcCall.options = new LinkedTreeMap<>();
        LinkedTreeMap<String, LinkedTreeMap<String, Object>> map = new LinkedTreeMap<>();
        funcCall.options.put("patterns", map);
    }

    public FuncCallBuilder addQuery(String query) {
        funcCall.func.query = query;
        return this;
    }

    public FuncCallBuilder addConnection(DataConnection connection) {
        funcCall.func.connection = connection;
        return this;
    }

    public FuncCallBuilder addFuncParams(List<FuncParam> funcParams) {
        funcCall.func.params.addAll(funcParams);
        return this;
    }

    public FuncCallBuilder addFuncParam(String type, String subType, String name, Object value, String patternType) {
        FuncParam funcParam = new FuncParam(type, name, value);
        funcParam.propertySubType = subType;
        funcParam.options = new LinkedTreeMap<>();
        funcParam.options.put("pattern", patternType);
        funcParam.options.put("default", "");
        funcCall.func.params.add(funcParam);
        return this;
    }

    @SuppressWarnings("unchecked")
    public <T> FuncCallBuilder addFuncCallOptionsPattern(String columnName, String expression,
                                                         String operator, Boolean include1,
                                                         Boolean include2, T... values) {
        LinkedTreeMap<String, Object> map1 = new LinkedTreeMap<>();
        map1.put("expression", expression);
        map1.put("colName", "");
        map1.put("values", Arrays.stream(values).filter(Objects::nonNull).map(Object::toString).collect(Collectors.toList()));
        map1.put("op", operator);
        map1.put("include1", include1);
        map1.put("include2", include2);
        LinkedTreeMap<String, LinkedTreeMap<String, Object>> patterns =
                (LinkedTreeMap<String, LinkedTreeMap<String, Object>>) funcCall.options.get("patterns");
        patterns.put(columnName, map1);
        return this;
    }

    public FuncCall build() {
        FuncCall returnFuncCall = funcCall;
        restore();
        return returnFuncCall;
    }
}

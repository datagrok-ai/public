package grok_connect.providers.utils;

import com.google.gson.internal.LinkedTreeMap;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataQuery;
import grok_connect.connectors_info.FuncCall;
import grok_connect.connectors_info.FuncParam;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

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

    public void restore() {
        funcCall = new FuncCall();
        funcCall.func = new DataQuery();
        funcCall.func.params = new ArrayList<>();
        funcCall.options = new LinkedTreeMap<String, Object>();
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

    public FuncCallBuilder addFuncParam(String type, String name, Object value, String patternType) {
        FuncParam funcParam = new FuncParam(type, name, value);
        funcParam.options = new LinkedTreeMap<>();
        funcParam.options.put("pattern", patternType);
        funcParam.options.put("default", "");
        funcCall.func.params.add(funcParam);
        return this;
    }

    @SuppressWarnings("unchecked")
    public <T> FuncCallBuilder addFuncCallOptionsPattern(String columnName, String expression,
                                                         String operator, T... values) {
        LinkedTreeMap<String, Object> map1 = new LinkedTreeMap<>();
        map1.put("expression", expression);
        map1.put("colName", "");
        map1.put("values", Arrays.stream(values).collect(Collectors.toList()));
        map1.put("op", operator);
        LinkedTreeMap<String, LinkedTreeMap<String, Object>> map2 = new LinkedTreeMap<>();
        map2.put(columnName, map1);
        funcCall.options.put("patterns", map2);
        return this;
    }

    public FuncCall build() {
        FuncCall returnFuncCall = funcCall;
        restore();
        return returnFuncCall;
    }
}

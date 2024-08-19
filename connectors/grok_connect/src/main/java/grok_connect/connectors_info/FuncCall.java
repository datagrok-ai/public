package grok_connect.connectors_info;

import java.util.*;


public class FuncCall {
    public static final String DEBUG_QUERY_KEY = "debug";
    public String id;
    public DataQuery func;
    public Map<String, Object> options;
    public Map<String, Object> parameterValues = new HashMap<>();
    public Map<String, Object> aux = new HashMap<>();
    public String log;
    public boolean debugQuery;

    public void setParamValues() {
        for (String paramName: parameterValues.keySet()) {
            for (FuncParam param: func.getInputParams()) {
                if (param.name.equals(paramName)) {
                    param.setValue(parameterValues.get(paramName));
                    break;
                }
            }
        }
    }

    public void afterDeserialization() {
        this.debugQuery = options.get(DEBUG_QUERY_KEY) != null && options.get(DEBUG_QUERY_KEY).equals(Boolean.TRUE);
    }
}

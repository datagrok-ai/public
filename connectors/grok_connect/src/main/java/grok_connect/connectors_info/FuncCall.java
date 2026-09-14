package grok_connect.connectors_info;

import java.util.*;


public class FuncCall {
    public static final String DEBUG_QUERY_KEY = "debug";
    public static final String LOG_QUERY_TEXT_KEY = "logQueryText";
    public String id;
    public DataQuery func;
    public Map<String, Object> options;
    public Map<String, Object> parameterValues = new HashMap<>();
    public Map<String, Object> aux = new HashMap<>();
    public String log;
    public boolean debugQuery;
    public boolean logQueryText;

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
        this.debugQuery = Boolean.TRUE.equals(options.get(DEBUG_QUERY_KEY));
        this.logQueryText = debugQuery || Boolean.TRUE.equals(options.get(LOG_QUERY_TEXT_KEY));
    }
}

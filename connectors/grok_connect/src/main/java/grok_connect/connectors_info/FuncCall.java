package grok_connect.connectors_info;

import java.util.*;


public class FuncCall
{
    public DataQuery func;
    public Map options;
    public Map<String, Object> parameterValues;

    public void setParamValues() {
        for (String paramName: parameterValues.keySet()) {
            for (FuncParam param: func.getInputParams()) {
                if (param.name.equals(paramName)) {
                    param.value = parameterValues.get(paramName);
                    break;
                }
            }
        }
    }
}

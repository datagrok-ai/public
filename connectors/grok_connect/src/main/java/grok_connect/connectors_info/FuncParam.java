package grok_connect.connectors_info;

import java.util.*;


public class FuncParam {
    public String propertyType;
    public String name;
    public Object value;
    public boolean isInput;
    public Map<String, String> options;

    public FuncParam() {
    }

    public FuncParam(String type, String name, Object value) {
        isInput = true;
        this.name = name;

        propertyType = type;
        if (type.equals("int") && (value instanceof Double)) {
            this.value = ((Double)value).intValue();
        } else {
            this.value = value;
        }
    }
}

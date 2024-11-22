package grok_connect.connectors_info;

import java.util.*;

import com.google.gson.annotations.SerializedName;
import serialization.Types;


public class FuncParam {
    @SerializedName("#type")
    public String type = "FuncParam";
    public String propertyType;
    public String propertySubType;
    public String name;
    public Object value;
    public boolean isInput = true;
    public Map<String, String> options;

    public FuncParam() {
    }

    public FuncParam(String type, String name, Object value) {
        this.name = name;
        propertyType = type;
        setValue(value);
    }

    public void setValue(Object value) {
        if (value instanceof String && !propertyType.equals(Types.STRING)) {
            if (propertyType.equals(Types.INT))
                this.value = Integer.parseInt((String)value);
            else if (propertyType.equals(Types.FLOAT))
                this.value = Float.parseFloat((String)value);
            else if (propertyType.equals(Types.BOOL))
                this.value = value.equals("true");
            else
                this.value = value;
        } else {
            if (propertyType.equals(Types.INT) && (value instanceof Double))
                this.value = ((Double) value).intValue();
            else
                this.value = value;
        }
    }

    @Override
    public String toString() {
        return "FuncParam{" +
                "propertyType='" + propertyType + '\'' +
                ", propertySubType='" + propertySubType + '\'' +
                ", name='" + name + '\'' +
                ", value=" + value +
                ", isInput=" + isInput +
                ", options=" + options +
                '}';
    }
}

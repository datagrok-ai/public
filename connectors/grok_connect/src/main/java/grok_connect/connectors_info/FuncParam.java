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
            switch (propertyType) {
                case Types.INT:
                    this.value = ((Double) Double.parseDouble((String) value)).intValue();
                    break;
                case Types.FLOAT:
                    this.value = Float.parseFloat((String) value);
                    break;
                case Types.BOOL:
                    this.value = value.equals("true");
                    break;
                default:
                    this.value = value;
                    break;
            }
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

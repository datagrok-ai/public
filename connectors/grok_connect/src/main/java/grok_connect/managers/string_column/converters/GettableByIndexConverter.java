package grok_connect.managers.string_column.converters;

import com.datastax.oss.driver.api.core.data.GettableByIndex;
import com.datastax.oss.driver.internal.core.data.DefaultUdtValue;
import grok_connect.managers.Converter;

public class GettableByIndexConverter implements Converter<String> {
    @Override
    public String convert(Object value) {
        return value instanceof GettableByIndex ? convertCassandraComplexType(value) :
                value.toString();
    }

    private String convertCassandraComplexType(Object value) {
        GettableByIndex complex = (GettableByIndex) value;
        StringBuilder builder = new StringBuilder("{");
        for (int i = 0; i < complex.size(); i++) {
            Object object = complex.getObject(i);
            if (object instanceof GettableByIndex) {
                builder.append(convertCassandraComplexType(object));
            } else {
                if (value instanceof DefaultUdtValue) {
                    String fieldName = ((DefaultUdtValue) value).getType()
                            .getFieldNames().get(i).toString();
                    builder.append(fieldName).append("=");
                }
                builder.append(object);
            }
            if (i != complex.size() - 1) {
                builder.append(", ");
            }
        }
        builder.append("}");
        return builder.toString();
    }
}

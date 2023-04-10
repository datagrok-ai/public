package grok_connect.converter.array.impl;

import grok_connect.converter.Converter;
import java.sql.Array;
import java.sql.SQLException;
import java.util.Arrays;

public class SQLArrayConverter implements Converter<String> {
    @Override
    public String convert(Object value) {
        if (value == null) {
            return Arrays.toString(new String[]{});
        }
        Array sqlArray = ((Array) value);
        try {
            Object[] array = (Object[]) sqlArray.getArray();
            return Arrays.toString(array);
        } catch (SQLException e) {
            throw new RuntimeException("Something went wrong when converting SQL ARRAY type", e);
        }
    }
}

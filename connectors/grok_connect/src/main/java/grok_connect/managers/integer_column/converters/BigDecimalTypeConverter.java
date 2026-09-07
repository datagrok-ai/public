package grok_connect.managers.integer_column.converters;

import grok_connect.managers.Converter;
import java.math.BigDecimal;

public class BigDecimalTypeConverter implements Converter<Integer> {
    @Override
    public Integer convert(Object value) {
        return ((BigDecimal) value).intValue();
    }
}

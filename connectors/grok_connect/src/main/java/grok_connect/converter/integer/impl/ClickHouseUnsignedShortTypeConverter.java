package grok_connect.converter.integer.impl;

import com.clickhouse.data.value.UnsignedShort;
import grok_connect.converter.Converter;

public class ClickHouseUnsignedShortTypeConverter implements Converter<Integer> {
    @Override
    public Integer convert(Object value) {
        return ((UnsignedShort) value).intValue();
    }
}

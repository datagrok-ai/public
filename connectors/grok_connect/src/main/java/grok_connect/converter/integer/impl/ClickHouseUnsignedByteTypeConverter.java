package grok_connect.converter.integer.impl;

import com.clickhouse.data.value.UnsignedByte;
import grok_connect.converter.Converter;

public class ClickHouseUnsignedByteTypeConverter implements Converter<Integer> {
    @Override
    public Integer convert(Object value) {
        return ((UnsignedByte) value).intValue();
    }
}

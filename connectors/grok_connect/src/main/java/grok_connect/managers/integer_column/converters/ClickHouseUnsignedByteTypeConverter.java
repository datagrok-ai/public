package grok_connect.managers.integer_column.converters;

import com.clickhouse.data.value.UnsignedByte;
import grok_connect.managers.Converter;

public class ClickHouseUnsignedByteTypeConverter implements Converter<Integer> {
    @Override
    public Integer convert(Object value) {
        return ((UnsignedByte) value).intValue();
    }
}

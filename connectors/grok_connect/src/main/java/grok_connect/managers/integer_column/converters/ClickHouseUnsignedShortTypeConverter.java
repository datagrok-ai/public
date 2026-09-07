package grok_connect.managers.integer_column.converters;

import com.clickhouse.data.value.UnsignedShort;
import grok_connect.managers.Converter;

public class ClickHouseUnsignedShortTypeConverter implements Converter<Integer> {
    @Override
    public Integer convert(Object value) {
        return ((UnsignedShort) value).intValue();
    }
}

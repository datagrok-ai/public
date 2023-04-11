package grok_connect.converter.integer.impl;

import com.clickhouse.data.value.UnsignedShort;
import grok_connect.converter.Converter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class ClickHouseUnsignedShortTypeConverter implements Converter<Integer> {
    private static final Logger LOGGER = LoggerFactory.getLogger(ClickHouseUnsignedShortTypeConverter.class);

    @Override
    public Integer convert(Object value) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        return ((UnsignedShort) value).intValue();
    }
}

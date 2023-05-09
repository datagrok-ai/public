package grok_connect.column.integer;

import com.clickhouse.data.value.UnsignedByte;
import com.clickhouse.data.value.UnsignedShort;
import grok_connect.column.ColumnManager;
import grok_connect.converter.Converter;
import grok_connect.converter.integer.BigDecimalTypeConverter;
import grok_connect.converter.integer.ByteTypeConverter;
import grok_connect.converter.integer.ClickHouseUnsignedByteTypeConverter;
import grok_connect.converter.integer.ClickHouseUnsignedShortTypeConverter;
import grok_connect.converter.integer.DoubleTypeConverter;
import grok_connect.converter.integer.FloatTypeConverter;
import grok_connect.converter.integer.LongTypeConverter;
import grok_connect.converter.integer.ShortTypeConverter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.Column;
import serialization.IntColumn;
import java.math.BigDecimal;
import java.util.HashMap;
import java.util.Map;

public class DefaultIntColumnManager implements ColumnManager<Integer> {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultIntColumnManager.class);
    private static final Converter<Integer> DEFAULT_CONVERTER = value -> (Integer) value;
    private final Map<Class<?>, Converter<Integer>> converterMap;

    {
        converterMap = new HashMap<>();
        converterMap.put(BigDecimal.class, new BigDecimalTypeConverter());
        converterMap.put(Byte.class, new ByteTypeConverter());
        converterMap.put(UnsignedByte.class, new ClickHouseUnsignedByteTypeConverter());
        converterMap.put(UnsignedShort.class, new ClickHouseUnsignedShortTypeConverter());
        converterMap.put(Double.class, new DoubleTypeConverter());
        converterMap.put(Float.class, new FloatTypeConverter());
        converterMap.put(Long.class, new LongTypeConverter());
        converterMap.put(Short.class, new ShortTypeConverter());
    }

    @Override
    public Integer convert(Object value, Object... args) {
        LOGGER.trace("convert method was called");
        if (value == null) {
            LOGGER.trace("value is null");
            return null;
        }
        Class<?> aClass = value.getClass();
        Converter<Integer> converter = converterMap
                .getOrDefault(aClass, DEFAULT_CONVERTER);
        return converter.convert(value);
    }

    @Override
    public boolean isApplicable(int type, String typeName, int precision, int scale) {
        return (type == java.sql.Types.INTEGER) || (type == java.sql.Types.TINYINT) || (type == java.sql.Types.SMALLINT)
                || typeName.equalsIgnoreCase("int4") || typeName.equalsIgnoreCase("int2")
                || typeName.equalsIgnoreCase("int") || typeName.equalsIgnoreCase("serial2")
                || typeName.equalsIgnoreCase("serial4") || typeName.equalsIgnoreCase("UInt16")
                || typeName.equalsIgnoreCase("UInt8") || (typeName.equalsIgnoreCase("NUMBER")
                && precision < 10 && scale == 0);
    }

    @Override
    public boolean isApplicable(Object o) {
        return o instanceof Byte || o instanceof Short || o instanceof Integer;
    }

    @Override
    public Column getColumn() {
        return new IntColumn();
    }

    @Override
    public Column getColumnWithInitSize(int size) {
        return new IntColumn(new Integer[size]);
    }
}

package grok_connect.managers.string_column.converters;

import grok_connect.managers.Converter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.sql.Array;
import java.sql.SQLException;
import java.util.Arrays;

public class SQLArrayTypeConverter implements Converter<String> {
    private static final Logger LOGGER = LoggerFactory.getLogger(SQLArrayTypeConverter.class);

    @Override
    public String convert(Object value) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        Array sqlArray = ((Array) value);
        try {
            Object[] array = (Object[]) sqlArray.getArray();
            return Arrays.toString(array);
        } catch (SQLException e) {
            LOGGER.debug("couldn't convert", e);
            throw new RuntimeException("Something went wrong when converting SQL ARRAY type", e);
        }
    }
}

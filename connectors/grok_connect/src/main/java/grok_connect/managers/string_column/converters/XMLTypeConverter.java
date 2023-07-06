package grok_connect.managers.string_column.converters;

import grok_connect.managers.Converter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.sql.SQLException;
import java.sql.SQLXML;

public class XMLTypeConverter implements Converter<String> {
    private static final Logger LOGGER = LoggerFactory.getLogger(XMLTypeConverter.class);

    @Override
    public String convert(Object value) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        try {
            SQLXML sqlxml = (SQLXML)value;
            return sqlxml.getString();
        } catch (SQLException e) {
            LOGGER.debug("couldn't convert", e);
            throw new RuntimeException("Something went wrong when converting xml to string", e);
        }
    }
}

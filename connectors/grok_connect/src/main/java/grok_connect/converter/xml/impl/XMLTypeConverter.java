package grok_connect.converter.xml.impl;

import grok_connect.converter.Converter;
import grok_connect.converter.time.impl.ZonedDateTimeTypeConverter;
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

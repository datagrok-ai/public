package grok_connect.converter.string.impl;

import grok_connect.converter.Converter;
import org.apache.commons.io.IOUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.Reader;
import java.io.StringWriter;
import java.sql.Clob;
import java.sql.SQLException;

public class ClobTypeConverter implements Converter<String> {
    private static final Logger LOGGER = LoggerFactory.getLogger(ClobTypeConverter.class);

    @Override
    public String convert(Object value) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        try {
            Reader reader = ((Clob)value).getCharacterStream();
            StringWriter writer = new StringWriter();
            IOUtils.copy(reader, writer);
            return writer.toString();
        } catch (SQLException | IOException e) {
            LOGGER.debug("couldn't convert", e);
            throw new RuntimeException("Something went wrong when converting clob type to string", e);
        }
    }
}

package grok_connect.managers.string_column.converters;

import grok_connect.managers.Converter;
import org.apache.commons.io.IOUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;
import java.sql.Blob;
import java.sql.SQLException;

public class BlobTypeConverter implements Converter<String> {
    private static final Logger LOGGER = LoggerFactory.getLogger(BlobTypeConverter.class);

    @Override
    public String convert(Object value) {
        if (value == null) {
            LOGGER.trace("value is null");
            return "";
        }
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        try {
            InputStream reader = ((Blob) value).getBinaryStream();
            StringWriter writer = new StringWriter();
            IOUtils.copy(reader, writer);
            return writer.toString();
        } catch (SQLException | IOException e) {
            LOGGER.debug("couldn't convert", e);
            throw new RuntimeException("Something went wrong when converting clob type to string", e);
        }
    }
}

package grok_connect.converter.string.impl;

import grok_connect.converter.Converter;
import org.apache.commons.io.IOUtils;
import java.io.IOException;
import java.io.Reader;
import java.io.StringWriter;
import java.sql.Clob;
import java.sql.SQLException;

public class ClobTypeConverter implements Converter<String> {
    @Override
    public String convert(Object value) {
        try {
            Reader reader = ((Clob)value).getCharacterStream();
            StringWriter writer = new StringWriter();
            IOUtils.copy(reader, writer);
            return writer.toString();
        } catch (SQLException | IOException e) {
            throw new RuntimeException("Something went wrong when converting clob type to string", e);
        }
    }
}

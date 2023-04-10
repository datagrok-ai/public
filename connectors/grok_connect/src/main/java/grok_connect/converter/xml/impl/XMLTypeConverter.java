package grok_connect.converter.xml.impl;

import grok_connect.converter.Converter;
import java.sql.SQLException;
import java.sql.SQLXML;

public class XMLTypeConverter implements Converter<String> {
    @Override
    public String convert(Object value) {
        try {
            SQLXML sqlxml = (SQLXML)value;
            return sqlxml.getString();
        } catch (SQLException e) {
            throw new RuntimeException("Something went wrong when converting xml to string", e);
        }
    }
}

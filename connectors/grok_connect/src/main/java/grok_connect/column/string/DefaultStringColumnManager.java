package grok_connect.column.string;

import com.ibm.db2.jcc.am.c9;
import com.ibm.db2.jcc.am.db;
import grok_connect.column.ColumnManager;
import grok_connect.converter.Converter;
import grok_connect.converter.array.ArrayTypeConverter;
import grok_connect.converter.array.SQLArrayConverter;
import grok_connect.converter.string.ClobTypeConverter;
import grok_connect.converter.xml.XMLTypeConverter;
import grok_connect.resultset.ColumnMeta;
import oracle.sql.ARRAY;
import oracle.xdb.XMLType;
import org.postgresql.jdbc.PgSQLXML;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.Column;
import serialization.StringColumn;
import java.sql.Array;
import java.sql.Clob;
import java.sql.SQLXML;
import java.util.HashMap;
import java.util.Map;

public class DefaultStringColumnManager implements ColumnManager<String> {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultStringColumnManager.class);
    private static final Converter<String> DEFAULT_CONVERTER = value -> value == null ? ""
            : value.toString();
    private final Map<Class<?>, Converter<String>> converterMap;

    {
        converterMap = new HashMap<>();
        Converter<String> clobConverter = new ClobTypeConverter();
        converterMap.put(Clob.class, clobConverter);
        converterMap.put(c9.class, clobConverter);
        Converter<String> xmlTypeConverter = new XMLTypeConverter();
        converterMap.put(SQLXML.class, xmlTypeConverter);
        converterMap.put(PgSQLXML.class, xmlTypeConverter);
        converterMap.put(db.class, xmlTypeConverter);
        converterMap.put(XMLType.class, xmlTypeConverter);
        Converter<String> sqlArrayConverter = new SQLArrayConverter();
        converterMap.put(Object.class, new ArrayTypeConverter());
        converterMap.put(Array.class, sqlArrayConverter);
        converterMap.put(ARRAY.class, sqlArrayConverter);
    }

    @Override
    public String convert(Object value, String columnLabel) {
        LOGGER.trace("convert method was called");
        if (value == null) return "";
        Class<?> aClass = value.getClass();
        Converter<String> converter;
        if (aClass.isArray()) {
            converter = converterMap
                    .get(Object.class);

        } else {
            converter = converterMap
                    .getOrDefault(aClass, DEFAULT_CONVERTER);
        }
        return converter.convert(value);
    }

    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        int type = columnMeta.getType();
        String typeName = columnMeta.getTypeName();
        int precision = columnMeta.getPrecision();
        int scale = columnMeta.getScale();
        LOGGER.trace("isApplicable method was called with parameters: {}, {}, {}, {}", type,
                typeName, precision, scale);
        return type == java.sql.Types.ARRAY ||
                typeName.equalsIgnoreCase("ARRAY") ||
                ((type == java.sql.Types.VARCHAR)|| (type == java.sql.Types.CHAR) ||
                (type == java.sql.Types.LONGVARCHAR) || (type == java.sql.Types.CLOB) ||
                        (type == java.sql.Types.NCLOB) ||
                typeName.equalsIgnoreCase("varchar") ||
                typeName.equalsIgnoreCase("nvarchar") ||
                typeName.equalsIgnoreCase("nchar") ||
                typeName.equalsIgnoreCase("ntext")) &&
                !typeName.equalsIgnoreCase("uuid") &&
                !typeName.equalsIgnoreCase("set") ||
                (type == java.sql.Types.SQLXML ||
                        typeName.equalsIgnoreCase("xml"));
    }

    @Override
    public boolean isApplicable(Object o) {
        return o instanceof String;
    }

    @Override
    public Column getColumn() {
        return new StringColumn();
    }
}

package grok_connect.converter.xml;

import grok_connect.converter.AbstractConverterManager;
import grok_connect.converter.Converter;
import grok_connect.converter.xml.impl.XMLTypeConverter;
import grok_connect.type.TypeChecker;
import oracle.xdb.XMLType;
import org.postgresql.jdbc.PgSQLXML;
import com.ibm.db2.jcc.am.db;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.sql.SQLXML;
import java.util.HashMap;
import java.util.Map;

public class XmlConverterManager extends AbstractConverterManager<String> {
    private static final Logger LOGGER = LoggerFactory.getLogger(XmlConverterManager.class);
    private static final Converter<String> defaultConverter = Object::toString;
    private final Map<Class<?>, Converter<String>> converterMap;

    public XmlConverterManager(TypeChecker typeChecker) {
        super(typeChecker);
        converterMap = new HashMap<>();
        XMLTypeConverter xmlTypeConverter = new XMLTypeConverter();
        converterMap.put(SQLXML.class, xmlTypeConverter);
        converterMap.put(PgSQLXML.class, xmlTypeConverter);
        converterMap.put(db.class, xmlTypeConverter);
        converterMap.put(XMLType.class, xmlTypeConverter);
    }

    @Override
    public String convert(Object value, Object...args) {
        LOGGER.trace("convert method was called");
        if (value == null) {
            LOGGER.trace("value is null");
            return null;
        }
        Class<?> aClass = value.getClass();
        Converter<String> converter = converterMap
                .get(aClass);
        if (converter != null) {
            LOGGER.trace("using defined xmlTypeConverter for class {}", aClass);
            return converter.convert(value);
        }
        LOGGER.debug("couldn't find xmlTypeConverter, using default converter for class {}",
                aClass);
        return defaultConverter.convert(value);
    }
}

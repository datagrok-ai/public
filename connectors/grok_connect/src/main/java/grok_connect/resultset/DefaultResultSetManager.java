package grok_connect.resultset;

import grok_connect.column.BigIntColumnProvider;
import grok_connect.column.BoolColumnProvider;
import grok_connect.column.ColumnProvider;
import grok_connect.column.ComplexTypeColumnProvider;
import grok_connect.column.DateTimeColumnProvider;
import grok_connect.column.FloatColumnProvider;
import grok_connect.column.IntColumnProvider;
import grok_connect.column.StringColumnProvider;
import grok_connect.converter.ConverterManager;
import grok_connect.converter.array.ArrayConverterManager;
import grok_connect.converter.bigint.BigIntConverterManager;
import grok_connect.converter.bitstring.BitStringConverterManager;
import grok_connect.converter.bool.BoolTypeConverterManager;
import grok_connect.converter.complex.ComplexTypeConverterManager;
import grok_connect.converter.float_type.FloatTypeConverterManager;
import grok_connect.converter.integer.IntegerTypeConverterManager;
import grok_connect.converter.string.StringTypeConverterManager;
import grok_connect.converter.time.TimeTypeConverterManager;
import grok_connect.converter.xml.XmlConverterManager;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.Column;
import serialization.StringColumn;
import java.util.ArrayList;
import java.util.List;

public class DefaultResultSetManager implements ResultSetManager {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultResultSetManager.class);
    private static final int COLUMN_NAME_INDEX = 0;
    private static ResultSetManager defaultManager;
    private final List<ConverterManager<?>> converterManagers;
    private final List<ColumnProvider> columnProviders;

    public DefaultResultSetManager(List<ConverterManager<?>> converterManagers,
                                   List<ColumnProvider> columnProviders) {
        this.columnProviders = columnProviders;
        this.converterManagers = converterManagers;
    }

    public static synchronized ResultSetManager getDefaultManager() {
        if (defaultManager == null) {
            defaultManager = new DefaultResultSetManager(getDefaultConverterManagers(),
                    getDefaultColumnProviders());
        }
        return defaultManager;
    }

    @SuppressWarnings("unchecked")
    @Override
    public <T> T convert(Object o, int type, String typeName, int precision, int scale, Object...args) {
        LOGGER.trace("convert method was called");
        for (ConverterManager<?> manager : converterManagers) {
            if (manager.isSupported(type, typeName, precision, scale)) {
                LOGGER.trace("found suitable converter manager");
                return (T) manager.convert(o, args.length == 0 ? null : args[COLUMN_NAME_INDEX]);
            }
        }
        LOGGER.trace("can't find suitable converter manager, return as a string");
        return o == null ? null : (T) o.toString();
    }

    @Override
    public Column getColumn(int type, String typeName, int precision, int scale) {
        LOGGER.trace("getColumn method was called");
        for (ColumnProvider columnProvider: columnProviders) {
            if (columnProvider.isSupported(type, typeName, precision, scale)) {
                LOGGER.trace("found suitable column provider");
                return columnProvider.get();
            }
        }
        LOGGER.trace("couldn't find suitable column, return StringColumn");
        return new StringColumn();
    }

    @Override
    public Column getColumnWithInitSize(int type, String typeName, int precision, int scale, int size) {
        LOGGER.trace("getColumn method with init size was called");
        for (ColumnProvider columnProvider: columnProviders) {
            if (columnProvider.isSupported(type, typeName, precision, scale)) {
                LOGGER.trace("found suitable column provider");
                return columnProvider.getWithInitSize(size);
            }
        }
        LOGGER.trace("couldn't find suitable column, return StringColumn");
        return new StringColumn();
    }

    @Override
    public Column getColumn(Object o) {
        LOGGER.trace("getColumn method with init size was called for object");
        for (ColumnProvider columnProvider: columnProviders) {
            if (columnProvider.isSupported(o)) {
                LOGGER.trace("found suitable column provider");
                return columnProvider.get();
            }
        }
        LOGGER.trace("couldn't find suitable column, return StringColumn");
        return new StringColumn();
    }

    public static List<ColumnProvider> getDefaultColumnProviders() {
        List<ColumnProvider> columnProviders = new ArrayList<>();
        columnProviders.add(new IntColumnProvider());
        columnProviders.add(new ComplexTypeColumnProvider());
        columnProviders.add(new BigIntColumnProvider());
        columnProviders.add(new StringColumnProvider());
        columnProviders.add(new FloatColumnProvider());
        columnProviders.add(new DateTimeColumnProvider());
        columnProviders.add(new BoolColumnProvider());
        return columnProviders;
    }

    public static List<ConverterManager<?>> getDefaultConverterManagers() {
        List<ConverterManager<?>> converterManagers = new ArrayList<>();
        converterManagers.add(new BigIntConverterManager());
        converterManagers.add(new IntegerTypeConverterManager());
        converterManagers.add(new ComplexTypeConverterManager());
        converterManagers.add(new ArrayConverterManager());
        converterManagers.add(new BitStringConverterManager());
        converterManagers.add(new BoolTypeConverterManager());
        converterManagers.add(new FloatTypeConverterManager());
        converterManagers.add(new StringTypeConverterManager());
        converterManagers.add(new TimeTypeConverterManager());
        converterManagers.add(new XmlConverterManager());
        return converterManagers;
    }
}

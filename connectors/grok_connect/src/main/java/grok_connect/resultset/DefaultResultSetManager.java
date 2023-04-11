package grok_connect.resultset;

import grok_connect.column.ColumnProvider;
import grok_connect.converter.ConverterManager;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.Column;
import serialization.StringColumn;
import java.util.List;

public class DefaultResultSetManager implements ResultSetManager {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultResultSetManager.class);
    private static final int COLUMN_NAME_INDEX = 0;
    private final List<ConverterManager<?>> converterManagers;
    private final List<ColumnProvider> columnProviders;

    public DefaultResultSetManager(List<ConverterManager<?>> converterManagers,
                                   List<ColumnProvider> columnProviders) {
        this.columnProviders = columnProviders;
        this.converterManagers = converterManagers;
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
}

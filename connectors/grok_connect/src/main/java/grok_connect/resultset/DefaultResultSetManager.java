package grok_connect.resultset;

import grok_connect.column.ColumnManager;
import grok_connect.column.bigint.DefaultBigIntColumnManager;
import grok_connect.column.bool.DefaultBoolColumnManager;
import grok_connect.column.complex.DefaultComplexColumnManager;
import grok_connect.column.datetime.DefaultDateTimeColumnManager;
import grok_connect.column.floats.DefaultFloatColumnManager;
import grok_connect.column.integer.DefaultIntColumnManager;
import grok_connect.column.string.DefaultStringColumnManager;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.Column;
import serialization.StringColumn;
import serialization.Types;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

public class DefaultResultSetManager implements ResultSetManager {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultResultSetManager.class);
    private static final int COLUMN_NAME_INDEX = 0;
    private static ResultSetManager defaultManager;
    private final Collection<ColumnManager<?>> columnManagers;

    public DefaultResultSetManager(Collection<ColumnManager<?>> columnManagers) {
        this.columnManagers = columnManagers;
    }

    public static synchronized ResultSetManager getDefaultManager() {
        if (defaultManager == null) {
            defaultManager = new DefaultResultSetManager(getDefaultManagersMap().values());
        }
        return defaultManager;
    }

    public static ResultSetManager fromManagersMap(Map<String, ColumnManager<?>> managers) {
        return new DefaultResultSetManager(managers.values());
    }

    public static Map<String, ColumnManager<?>> getDefaultManagersMap() {
        Map<String, ColumnManager<?>> map = new HashMap<>();
        map.put(Types.BIG_INT, new DefaultBigIntColumnManager());
        map.put(Types.BOOL, new DefaultBoolColumnManager());
        map.put(Types.COLUMN_LIST, new DefaultComplexColumnManager());
        map.put(Types.DATE_TIME, new DefaultDateTimeColumnManager());
        map.put(Types.FLOAT, new DefaultFloatColumnManager());
        map.put(Types.INT, new DefaultIntColumnManager());
        map.put(Types.STRING, new DefaultStringColumnManager());
        return map;
    }

    @SuppressWarnings("unchecked")
    @Override
    public <T> T convert(Object o, int type, String typeName, int precision, int scale, Object...args) {
        LOGGER.trace("convert method was called");
        for (ColumnManager<?> manager : columnManagers) {
            if (manager.isApplicable(type, typeName, precision, scale)) {
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
        for (ColumnManager<?> manager : columnManagers) {
            if (manager.isApplicable(type, typeName, precision, scale)) {
                LOGGER.trace("found suitable column provider");
                return manager.getColumn();
            }
        }
        LOGGER.trace("couldn't find suitable column, return StringColumn");
        return new StringColumn();
    }

    @Override
    public Column getColumnWithInitSize(int type, String typeName, int precision, int scale, int size) {
        LOGGER.trace("getColumnWithInitSize method was called");
        for (ColumnManager<?> manager : columnManagers) {
            if (manager.isApplicable(type, typeName, precision, scale)) {
                LOGGER.trace("found suitable column provider");
                return manager.getColumnWithInitSize(size);
            }
        }
        LOGGER.trace("couldn't find suitable column, return StringColumn");
        return new StringColumn(new String[size]);
    }

    @Override
    public Column getColumn(Object o) {
        LOGGER.trace("getColumn method was called");
        for (ColumnManager<?> manager : columnManagers) {
            if (manager.isApplicable(o)) {
                LOGGER.trace("found suitable column provider");
                return manager.getColumn();
            }
        }
        LOGGER.trace("couldn't find suitable column, return StringColumn");
        return new StringColumn();
    }
}

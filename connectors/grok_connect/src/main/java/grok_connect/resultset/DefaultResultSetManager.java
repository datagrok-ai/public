package grok_connect.resultset;

import grok_connect.managers.ColumnManager;
import grok_connect.managers.bigint_column.DefaultBigIntColumnManager;
import grok_connect.managers.bool_column.DefaultBoolColumnManager;
import grok_connect.managers.complex_column.DefaultComplexColumnManager;
import grok_connect.managers.datetime_column.DefaultDateTimeColumnManager;
import grok_connect.managers.float_column.DefaultFloatColumnManager;
import grok_connect.managers.integer_column.DefaultIntColumnManager;
import grok_connect.managers.string_column.DefaultStringColumnManager;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.Column;
import serialization.StringColumn;
import serialization.Types;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class DefaultResultSetManager implements ResultSetManager {
    protected final Logger logger = LoggerFactory.getLogger(this.getClass().getName());
    protected final Collection<ColumnManager<?>> columnManagers;
    protected final List<Column> columns;
    protected final List<ColumnMeta> columnsMeta;
    protected final List<ColumnManager<?>> currentManagers;

    public DefaultResultSetManager(Collection<ColumnManager<?>> columnManagers) {
        this.columnManagers = columnManagers;
        columns = new ArrayList<>();
        columnsMeta = new ArrayList<>();
        currentManagers = new ArrayList<>();
    }

    public static ResultSetManager getDefaultManager() {
        return fromManagersMap(getDefaultManagersMap());
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
    public <T> T convert(Object o, ColumnMeta columnMeta) {
        logger.trace("convert method was called");
        for (ColumnManager<?> manager : columnManagers) {
            if (manager.isApplicable(columnMeta)) {
                logger.trace("found suitable converter manager");
                currentManagers.add(manager);
                return (T) manager.convert(o, columnMeta.getColumnLabel());
            }
        }
        logger.trace("can't find suitable converter manager, return as a string");
        return o == null ? null : (T) o.toString();
    }

    @Override
    public Column getColumn(ColumnMeta columnMeta) {
        logger.trace("getColumn method was called");
        for (ColumnManager<?> manager : columnManagers) {
            if (manager.isApplicable(columnMeta)) {
                logger.trace("found suitable column provider");
                return manager.getColumn();
            }
        }
        logger.trace("couldn't find suitable column, return StringColumn");
        return new StringColumn();
    }

    @Override
    public Column getColumn(Object o) {
        logger.trace("getColumn method was called");
        for (ColumnManager<?> manager : columnManagers) {
            if (manager.isApplicable(o)) {
                logger.trace("found suitable column provider");
                return manager.getColumn();
            }
        }
        logger.trace("couldn't find suitable column, return StringColumn");
        return new StringColumn();
    }

    @Override
    public void processValue(Object o, int index, ResultSetMetaData meta) {
        boolean inBounds = (index - 1 >= 0) && (index - 1 < columns.size());
        if (!inBounds) {
            processHeaders(o, index, meta);
        } else {
            columns.get(index - 1).add(currentManagers.get(index - 1)
                    .convert(o, columnsMeta.get(index - 1).getColumnLabel()));
        }
    }

    @Override
    public List<Column> getProcessedColumns() {
        return columns;
    }

    protected void processHeaders(Object o, int index, ResultSetMetaData meta) {
        ColumnMeta columnMeta = getColumnMeta(index, meta);
        columnsMeta.add(index - 1, columnMeta);
        Column column = getColumn(columnMeta);
        column.name = columnMeta.getColumnLabel();
        column.add(convert(o, columnMeta));
        columns.add(index - 1, column);
    }

    private ColumnMeta getColumnMeta(int index, ResultSetMetaData meta) {
        try {
            int type = meta.getColumnType(index);
            String typeName = meta.getColumnTypeName(index);
            int precision = meta.getPrecision(index);
            int scale = meta.getScale(index);
            String columnLabel = meta.getColumnLabel(index);
            return new ColumnMeta(type, typeName, precision, scale, columnLabel);
        } catch (SQLException e) {
            throw new RuntimeException("Something went wrong when processing ResultSetMetaData", e);
        }
    }
}

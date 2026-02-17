package grok_connect.resultset;

import grok_connect.log.EventType;
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
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

public class DefaultResultSetManager implements ResultSetManager {
    protected final Logger logger = LoggerFactory.getLogger(this.getClass().getName());
    protected final Collection<ColumnManager<?>> columnManagers;
    protected Column<?>[] columns;
    protected ColumnMeta[] columnsMeta;
    protected ColumnManager<?>[] currentManagers;
    protected boolean isInit = false;

    public DefaultResultSetManager(Collection<ColumnManager<?>> columnManagers) {
        this.columnManagers = columnManagers;
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

    @Override
    public void init(ResultSetMetaData meta, int initColumnSize) {
        try {
            int columnCount = meta.getColumnCount();
            columns = new Column<?>[columnCount];
            columnsMeta = new ColumnMeta[columnCount];
            currentManagers = new ColumnManager[columnCount];
            for (int i = 0; i < columnCount; i++) {
                ColumnMeta columnMeta = getColumnMeta(i + 1, meta, initColumnSize);
                columnsMeta[i] = columnMeta;
                ColumnManager<?> applicableColumnManager = getApplicableColumnManager(columnMeta);
                currentManagers[i] = applicableColumnManager;
                String columnName = columnMeta.getColumnLabel();
                columns[i] = applicableColumnManager.getColumn(columnName, initColumnSize);
            }
            isInit = true;
        } catch (SQLException e) {
            throw new RuntimeException("Something went wrong when init ResultSetManager", e);
        }
    }

    @Override
    public void processValue(Object o, int index) {
        if (isInit) {
            Column<?> current = columns[index - 1];
            try {
                setValue(o, index);
            } catch (RuntimeException e) {
                logger.error(EventType.ERROR.getMarker(),
                        "Couldn't convert object to be applicable for column {} of type {}. StringColumn fallback will be used",
                        columnsMeta[index - 1].getColumnLabel(), columns[index - 1].getType(), e);
                StringColumn stringColumn = new StringColumn(current.getName(), columnsMeta[index - 1].getColumnSize());
                ColumnManager<String> newManager = new DefaultStringColumnManager();
                for (int i = 0; i < current.getLength(); i++) {
                    Object value = current.get(i);
                    stringColumn.add(newManager.convert(value, columnsMeta[index - 1]));
                }
                stringColumn.add(newManager.convert(o, columnsMeta[index - 1]));
                columns[index - 1] = stringColumn;
                currentManagers[index - 1] = newManager;
            }
        }
        else
            throw new RuntimeException("ResultSetManager should be init");
    }

    @Override
    public Column<?>[] getProcessedColumns() {
        return columns;
    }

    @Override
    public void empty(int newColSize) {
        for (Column<?> col : columns) {
            col.setInitColumnSize(newColSize);
            col.empty();
        }
    }

    @Override
    public ColumnManager<?> getApplicableColumnManager(ColumnMeta meta) {
        for (ColumnManager<?> manager : columnManagers) {
            if (manager.isApplicable(meta)) {
                return manager;
            }
        }
        return new DefaultStringColumnManager();
    }

    @SuppressWarnings("unchecked")
    protected void setValue(Object o, int index) {
        ((Column<Object>) columns[index - 1]).add(currentManagers[index - 1].convert(o, columnsMeta[index - 1]));
    }

    protected ColumnMeta getColumnMeta(int index, ResultSetMetaData meta, int initColumnSize) {
        try {
            int type = meta.getColumnType(index);
            String typeName = meta.getColumnTypeName(index);
            int precision = meta.getPrecision(index);
            int scale = meta.getScale(index);
            String columnLabel = meta.getColumnLabel(index);
            return new ColumnMeta(type, typeName, precision, scale, columnLabel, initColumnSize);
        } catch (SQLException e) {
            throw new RuntimeException("Something went wrong when processing ResultSetMetaData", e);
        }
    }
}

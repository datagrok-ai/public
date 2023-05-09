package grok_connect.resultset;

import serialization.Column;
import java.sql.ResultSetMetaData;
import java.util.List;

public interface ResultSetManager {
    <T> T convert(Object o, ColumnMeta columnMeta);

    Column getColumn(ColumnMeta columnMeta);

    Column getColumn(Object o);

    void processValue(Object o, int index, ResultSetMetaData meta);

    List<Column> getProcessedColumns();
}

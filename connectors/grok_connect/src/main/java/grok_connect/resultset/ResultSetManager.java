package grok_connect.resultset;

import serialization.Column;
import java.sql.ResultSetMetaData;

public interface ResultSetManager {

    /**
     * Used to init all necessary managers, columns and columnsMetas in arrays for
     * fast access to them lately in processValue method.
     *
     * @param meta ResultSetMetaData (JDBC interface)
     */
    void init(ResultSetMetaData meta);

    /**
     * Used by ComplexTypeColumnManager. In future may be deleted.
     *
     * @param o Object to convert
     * @param columnMeta ColumnMeta for Object
     * @return converted Object
     */
    Object convert(Object o, ColumnMeta columnMeta);

    /**
     * Used by ComplexTypeColumnManager. In future may be deleted.
     *
     * @param o Object for which provide Column
     * @return Column for object
     */
    Column getColumn(Object o);

    /**
     * Used to process Object from ResultSet based on previously generated ColumnsMeta, Column and ColumnManager.
     *
     * @param o Object from ResultSet
     * @param index Column index
     */
    void processValue(Object o, int index);

    /**
     * Return filled Columns
     *
     * @return Column array
     */
    Column[] getProcessedColumns();

    /**
     * Empty all columns for reuse without the need to init ResultSetManager again
     */
    void empty();
}

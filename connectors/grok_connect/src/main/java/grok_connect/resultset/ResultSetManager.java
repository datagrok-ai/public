package grok_connect.resultset;

import grok_connect.managers.ColumnManager;
import serialization.Column;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;

public interface ResultSetManager {

    /**
     * Used to init all necessary managers, columns and columnsMetas in arrays for
     * fast access to them lately in processValue method.
     *
     * @param meta ResultSetMetaData (JDBC interface)
     */
    void init(ResultSetMetaData meta, int initColumnSize);

    /** Lets BigIntColumn report itself as int while its values fit; set before init. */
    void setInt8AsInt32(boolean value);

    boolean readFast(ResultSet resultSet, int index) throws SQLException;

    ColumnManager<?> getApplicableColumnManager(ColumnMeta meta);

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
    Column<?>[] getProcessedColumns();

    /**
     * Replaces the filled columns with fresh empty ones of the same classes without the need to init
     * ResultSetManager again; the previous columns stay intact for the DataFrame that holds them.
     */
    void detach(int nextColSize);
}

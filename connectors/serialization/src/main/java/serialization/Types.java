package serialization;

import java.util.*;

public class Types {
    public static final String INT = "int";
    public static final String FLOAT = "double";
    public static final String BOOL = "bool";
    public static final String STRING = "string";
    public static final String DATE_TIME = "datetime";
    public static final String BIG_INT = "bigint";
    public static final String NUM = "num";
    public static final String STRING_LIST = "string_list";
    public static final String OBJECT = "object";
    public static final String DATA_FRAME = "dataframe";
    public static final String DATA_FRAME_LIST = "dataframe_list";
    public static final String CELL = "cell";
    public static final String COLUMN = "column";
    public static final String COLUMN_LIST = "column_list";
    public static final String GRAPHICS = "graphics";
    public static final String ROW_FILTER = "tablerowfiltercall";
    public static final String COLUMN_FILTER = "colfiltercall";
    public static final String BIT_SET = "bitset";
    public static final String MAP = "map";
    public static final String DYNAMIC = "dynamic";
    public static final String LIST = "list";
    public static final String BLOB = "blob";

    /// Types supported by [DataFrame]
    public static final List<String> dataFrameColumnTypes = new ArrayList<String>() {{
        add(INT);
        add(FLOAT);
        add(BOOL);
        add(STRING);
        add(DATE_TIME);
        add(OBJECT);
    }};

    /// Numeric types supported by [DataFrame]
    public static final List<String> dataFrameNumericTypes = new ArrayList<String>() {{
        add(INT);
        add(FLOAT);
    }};
}

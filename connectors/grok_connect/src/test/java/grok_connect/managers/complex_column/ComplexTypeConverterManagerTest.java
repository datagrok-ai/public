package grok_connect.managers.complex_column;

import grok_connect.GrokConnect;
import grok_connect.managers.ColumnManager;
import grok_connect.resultset.ColumnMeta;
import grok_connect.utils.ProviderManager;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Assertions;
import serialization.Column;
import serialization.ComplexTypeColumn;
import serialization.DataFrame;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

class ComplexColumnManagerTest {
    private static final String firstColumnName = "id";
    private static final String secondColumnName = "firstName";
    private static final String thirdColumnName = "someList";
    private static final String fourthColumnName = "some_double";
    private static final String fifthColumnName = "integer";
    private static final String sixthColumnName = "bool";
    private static final String nestedColumnName = "nested";
    private static final Long firstColumnData = 1L;
    private static final String secondColumnData = "John";
    private static final List<String> thirdColumnData = new ArrayList<>();
    private static final Double fourthColumnData = 3.14;
    private static final Integer fifthColumnData = 2023;
    private static final Boolean sixthColumnData = true;
    private static ColumnManager<Map<String, Object>> complexColumnManager;
    private static ColumnMeta columnMeta;

    @BeforeAll
    public static void init() {
        thirdColumnData.add("work");
        thirdColumnData.add("rest");
        thirdColumnData.add("play");
        GrokConnect.providerManager = new ProviderManager();
        complexColumnManager = new DefaultComplexColumnManager();
        columnMeta = new ColumnMeta(1111, "test", -1, -1, "Column", 100);
    }

    @Test
    public void testNullArgument_ok() {
        Map<String, Object> result = Assertions.assertDoesNotThrow(() -> complexColumnManager
                .convert(null, columnMeta));
        Assertions.assertNotNull(result);
        Assertions.assertTrue(result.isEmpty());
    }

    @Test
    public void testFlatMapArgument_ok() {
        Map<String, Object> map = getFlatMap();
        Map<String, Object> converted = complexColumnManager.convert(map, columnMeta);
        Assertions.assertNotNull(converted);
        Assertions.assertEquals(map.size(), converted.size());
    }

    @Test
    public void testComplexTypeColumn_flatMap() {
        ComplexTypeColumn ctc = new ComplexTypeColumn("Column");
        ctc.add(getFlatMap());
        DataFrame df = new DataFrame();
        df.addColumns(new Column<?>[]{ctc});
        // CTC is expanded into flat columns at add time
        Assertions.assertTrue(df.getColumnCount() > 0);
        // Verify column names are prefixed with ComplexTypeColumn name
        for (Column<?> col : df.getColumns())
            Assertions.assertTrue(col.getName().startsWith("Column."), "Expected prefix 'Column.' but got: " + col.getName());
    }

    @Test
    public void testComplexTypeColumn_nestedMap() {
        ComplexTypeColumn ctc = new ComplexTypeColumn("Column");
        ctc.add(getNestedMap());
        DataFrame df = new DataFrame();
        df.addColumns(new Column<?>[]{ctc});
        // Should have flattened the nested map with dot-notation
        List<String> names = new ArrayList<>();
        for (Column<?> col : df.getColumns())
            names.add(col.getName());
        boolean hasNestedKey = false;
        for (String colName : names)
            if (colName.contains(nestedColumnName + ".")) {
                hasNestedKey = true;
                break;
            }
        Assertions.assertTrue(hasNestedKey, "Expected nested keys with dot notation, got: " + names);
    }

    @Test
    public void testComplexTypeColumn_multipleRows() {
        ComplexTypeColumn ctc = new ComplexTypeColumn("Column");
        ctc.add(getFlatMap());
        ctc.add(getFlatMap());
        Assertions.assertEquals(2, ctc.getLength());
        DataFrame df = new DataFrame();
        df.addColumns(new Column<?>[]{ctc});
        Assertions.assertTrue(df.getColumnCount() > 0);
        // Each column should have 2 rows
        for (Column<?> col : df.getColumns())
            Assertions.assertEquals(2, col.getLength(), "Column " + col.getName() + " should have 2 rows");
    }

    private Map<String, Object> getFlatMap() {
        Map<String, Object> argument = new HashMap<>();
        argument.put(firstColumnName, firstColumnData);
        argument.put(secondColumnName, secondColumnData);
        argument.put(thirdColumnName, thirdColumnData);
        argument.put(fourthColumnName, fourthColumnData);
        argument.put(fifthColumnName, fifthColumnData);
        argument.put(sixthColumnName, sixthColumnData);
        return argument;
    }

    private Map<String, Object> getNestedMap() {
        Map<String, Object> flatMap = getFlatMap();
        flatMap.put(nestedColumnName, getFlatMap());
        return flatMap;
    }
}

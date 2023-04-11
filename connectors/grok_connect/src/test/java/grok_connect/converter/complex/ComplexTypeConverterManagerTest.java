package grok_connect.converter.complex;

import grok_connect.GrokConnect;
import grok_connect.providers.utils.DataFrameComparator;
import grok_connect.type.complex.DefaultComplexTypeChecker;
import grok_connect.utils.ProviderManager;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Assertions;
import serialization.BigIntColumn;
import serialization.BoolColumn;
import serialization.Column;
import serialization.DateTimeColumn;
import serialization.FloatColumn;
import serialization.IntColumn;
import serialization.StringColumn;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

class ComplexTypeConverterManagerTest {
    private static final String DEFAULT_MAIN_COLUMN_NAME = "Column";
    private static final String firstColumnName = "id";
    private static final String secondColumnName = "firstName";
    private static final String thirdColumnName = "someList";
    private static final String fourthColumnName = "date";
    private static final String fifthColumnName = "some_double";
    private static final String sixthColumnName = "integer";
    private static final String seventhColumnName = "bool";
    private static final String nestedColumnName = "nested";
    private static final Long firstColumnData = 1L;
    private static final String secondColumnData = "John";
    private static final List<String> thirdColumnData = new ArrayList<>();
    private static final Date fourthColumnData = Date.from(Instant.now());
    private static final Double fifthColumnData = Double.NEGATIVE_INFINITY;
    private static final Integer sixthColumnData = 2023;
    private static final Boolean seventhColumnData = true;
    private static ComplexTypeConverterManager complexTypeConverterManager;
    private static DataFrameComparator dataFrameComparator;

    @BeforeAll
    public static void init() {
        dataFrameComparator = new DataFrameComparator();
        thirdColumnData.add("work");
        thirdColumnData.add("rest");
        thirdColumnData.add("play");
        ProviderManager providerManager = new ProviderManager();
        GrokConnect.providerManager = providerManager;
        complexTypeConverterManager = new ComplexTypeConverterManager(new DefaultComplexTypeChecker());
    }

    @Test
    public void testNullArgument_ok() {
        Assertions.assertDoesNotThrow(() -> complexTypeConverterManager
                .convert(null, DEFAULT_MAIN_COLUMN_NAME));
    }

    @Test
    public void testFlatMapArgument_ok() {
        List<Column> actual = Assertions.assertDoesNotThrow(() -> complexTypeConverterManager
                .convert(getFlatMap(), DEFAULT_MAIN_COLUMN_NAME));
        Assertions.assertTrue(dataFrameComparator.isColumnsEqualUnOrdered(getExpectedForFlatMap(), actual));
    }

    @Test
    public void testNestedMapArgument_ok() {
        List<Column> actual = Assertions.assertDoesNotThrow(() -> complexTypeConverterManager
                .convert(getNestedMap(), DEFAULT_MAIN_COLUMN_NAME));
        Assertions.assertTrue(dataFrameComparator.isColumnsEqualUnOrdered(getExpectedForNestedMap(), actual));
    }

    private Map<String, Object> getFlatMap() {
        Map<String, Object> argument = new HashMap<>();
        argument.put(firstColumnName, firstColumnData);
        argument.put(secondColumnName, secondColumnData);
        argument.put(thirdColumnName, thirdColumnData);
        argument.put(fourthColumnName, fourthColumnData);
        argument.put(fifthColumnName, fifthColumnData);
        argument.put(sixthColumnName, sixthColumnData);
        argument.put(seventhColumnName, seventhColumnData);
        return argument;
    }

    private Map<String, Object> getNestedMap() {
        Map<String, Object> flatMap = getFlatMap();
        flatMap.put(nestedColumnName, getFlatMap());
        return flatMap;
    }

    private Map<String, Object> getFlatMapWithNulls() {
        Map<String, Object> argument = new HashMap<>();
        argument.put(firstColumnName, null);
        argument.put(secondColumnName, null);
        argument.put(thirdColumnName, null);
        argument.put(fourthColumnName, null);
        argument.put(fifthColumnName, null);
        argument.put(sixthColumnName, null);
        argument.put(seventhColumnName, null);
        return argument;
    }

    private List<Column> getExpectedForFlatMap() {
        String formatName = "%s.%s";
        Column firstColumn = new BigIntColumn();
        firstColumn.add(firstColumnData.toString());
        firstColumn.name = String.format(formatName, DEFAULT_MAIN_COLUMN_NAME, firstColumnName);
        Column secondColumn = new StringColumn();
        secondColumn.add(secondColumnData);
        secondColumn.name = String.format(formatName, DEFAULT_MAIN_COLUMN_NAME, secondColumnName);
        Column thirdColumn = new StringColumn();
        thirdColumn.add(thirdColumnData.toString());
        thirdColumn.name = String.format(formatName, DEFAULT_MAIN_COLUMN_NAME, thirdColumnName);
        Column fourthColumn = new DateTimeColumn();
        fourthColumn.add(fourthColumnData.getTime() * 1000.0);
        fourthColumn.name = String.format(formatName, DEFAULT_MAIN_COLUMN_NAME, fourthColumnName);
        Column fifthColumn = new FloatColumn();
        fifthColumn.add(fifthColumnData.floatValue());
        fifthColumn.name = String.format(formatName, DEFAULT_MAIN_COLUMN_NAME, fifthColumnName);
        Column sixthColumn = new IntColumn();
        sixthColumn.add(sixthColumnData);
        sixthColumn.name = String.format(formatName, DEFAULT_MAIN_COLUMN_NAME, sixthColumnName);
        Column seventhColumn = new BoolColumn();
        seventhColumn.add(seventhColumnData);
        seventhColumn.name = String.format(formatName, DEFAULT_MAIN_COLUMN_NAME, seventhColumnName);
        List<Column> expected = new ArrayList<>();
        expected.addAll(Arrays.asList(firstColumn, secondColumn, thirdColumn, fourthColumn, fifthColumn,
                sixthColumn, seventhColumn));
        return expected;
    }

    private List<Column> getExpectedForNestedMap() {
        List<Column> expectedForFlatMap = getExpectedForFlatMap();
        List<Column> nested = getExpectedForFlatMap();
        nested.forEach(column -> {
            String[] split = column.name.split("\\.");
            column.name = String.format("%s.%s.%s", split[0], nestedColumnName, split[1]);

        });
        expectedForFlatMap.addAll(nested);
        return expectedForFlatMap;
    }
}

package grok_connect.providers.arguments_provider;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.DataFrameBuilder;
import grok_connect.providers.utils.DateParser;
import grok_connect.providers.utils.FuncCallBuilder;
import grok_connect.providers.utils.Parser;
import org.junit.jupiter.api.Named;
import org.junit.jupiter.params.provider.Arguments;
import serialization.BigIntColumn;
import serialization.DataFrame;
import serialization.DateTimeColumn;
import serialization.FloatColumn;
import serialization.IntColumn;
import serialization.StringColumn;
import java.math.BigInteger;
import java.nio.charset.StandardCharsets;
import java.util.stream.Stream;

@SuppressWarnings("unused")
public class SnowflakeObjectsMother {
    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new StringColumn(new String[] {"INFORMATION_SCHEMA", "PUBLIC"}),
                        "TABLE_SCHEMA")
                .build();
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> getSchema_ok() {
        String firstColumnName = "TABLE_SCHEMA";
        String secondColumnName = "TABLE_NAME";
        String thirdColumnName = "COLUMN_NAME";
        String fourthColumnName = "DATA_TYPE";
        String fifthColumnName = "IS_VIEW";
        String catalog = "TEST";
        String schema = "PUBLIC";
        String table = "MOCK_DATA";
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(10)
                .setColumn(new StringColumn(), "TABLE_CATALOG", new String[] {catalog, catalog,
                        catalog, catalog, catalog, catalog, catalog,
                        catalog, catalog, catalog})
                .setColumn(new StringColumn(), firstColumnName, new String[] {schema, schema,
                        schema, schema, schema, schema, schema,
                        schema, schema, schema})
                .setColumn(new StringColumn(), secondColumnName, new String[] {table, table,
                        table, table, table, table,
                        table, table, table, table})
                .setColumn(new StringColumn(), thirdColumnName, new String[] {"ID", "FIRST_NAME", "LAST_NAME", "EMAIL", "GENDER",
                        "IP_ADDRESS", "BOOL", "COUNTRY", "DATE", "SOME_NUMBER"})
                .setColumn(new StringColumn(), fourthColumnName, new String[] {"NUMBER", "TEXT",
                        "TEXT", "TEXT", "TEXT", "TEXT",
                        "BOOLEAN", "TEXT", "DATE", "FLOAT"})
                .setColumn(new IntColumn(), fifthColumnName, new Integer[] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0})
                .build();
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_dateTypes_ok() {
        Parser parser = new DateParser();
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd",
                        "2011-10-29")}), "D")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss.SS",
                        "2020-03-12 01:02:03.123")}), "T")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("HH:mm:ss.SSS X",
                        "04:05:06.789 +00:00")}), "TM")
                .build();
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM dates_table");
        return Stream.of(Arguments.of(Named.of("DATE TYPES SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_numericTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(4)
                .setColumn(new BigIntColumn(new String[]{"99999999999999999999999999999999999999",
                        "-99999999999999999999999999999999999999", "9", "10000"}), "NUM")
                .setColumn(new FloatColumn(new Float[]{255.5f, 15.5f, 255.2f, 44444.4f}), "NUM10")
                .setColumn(new FloatColumn(new Float[]{Float.NaN, Float.POSITIVE_INFINITY,
                        Float.NEGATIVE_INFINITY, 1.234E+2f}), "FL1")
                .setColumn(new FloatColumn(new Float[]{0.2f, 1.34f, 317f, 124.412f}), "FL2")
                .setColumn(new IntColumn(new Integer[]{100, 51200, 7780, 12}), "I")
                .build();
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM test_number");
        return Stream.of(Arguments.of(Named.of("NUMERIC TYPES SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_binaryType_ok() {
        String expectedHexRepr = String.format("%010X", new BigInteger(1, "Datagrok".getBytes(StandardCharsets.UTF_8)));
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[]{expectedHexRepr}), "B")
                .build();
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT to_char(B) AS B FROM binary_table");
        return Stream.of(Arguments.of(Named.of("BINARY TYPE SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_geoType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new BigIntColumn(new String[]{"1", "2"}), "ID")
                .setColumn(new StringColumn(new String[]{"{\n" +
                        "  \"coordinates\": [\n" +
                        "    -122.35,\n" +
                        "    37.55\n" +
                        "  ],\n" +
                        "  \"type\": \"Point\"\n" +
                        "}",
                        "{\n" +
                                "  \"coordinates\": [\n" +
                                "    [\n" +
                                "      -124.2,\n" +
                                "      42\n" +
                                "    ],\n" +
                                "    [\n" +
                                "      -120.01,\n" +
                                "      41.99\n" +
                                "    ]\n" +
                                "  ],\n" +
                                "  \"type\": \"LineString\"\n" +
                                "}"}), "G")
                .build();
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM geospatial_table");
        return Stream.of(Arguments.of(Named.of("GEO TYPE SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_semiStructuredTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new BigIntColumn(new String[]{"1", "2"}), "ID")
                .setColumn(new StringColumn(new String[]{"[\n" +
                        "  1,\n" +
                        "  2,\n" +
                        "  3\n" +
                        "]", "[\n" +
                        "  1,\n" +
                        "  2,\n" +
                        "  3,\n" +
                        "  undefined\n" +
                        "]"}), "ARRAY1")
                .setColumn(new StringColumn(new String[]{"{\n" +
                        "  \"key1\": \"value1\",\n" +
                        "  \"key2\": \"value2\"\n" +
                        "}",
                        "{\n" +
                                "  \"key1\": \"value1\",\n" +
                                "  \"key2\": null\n" +
                                "}"}), "VARIANT1")
                .setColumn(new StringColumn(new String[]{"{\n" +
                        "  \"outer_key1\": {\n" +
                        "    \"inner_key1A\": \"1a\",\n" +
                        "    \"inner_key1B\": \"1b\"\n" +
                        "  },\n" +
                        "  \"outer_key2\": {\n" +
                        "    \"inner_key2\": 2\n" +
                        "  }\n" +
                        "}", "{\n" +
                        "  \"outer_key1\": {\n" +
                        "    \"inner_key1A\": \"1a\",\n" +
                        "    \"inner_key1B\": null\n" +
                        "  },\n" +
                        "  \"outer_key2\": {\n" +
                        "    \"inner_key2\": 2\n" +
                        "  }\n" +
                        "}"}), "OBJECT1")
                .build();
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM demonstration1");
        return Stream.of(Arguments.of(Named.of("SEMI STRUCTURED TYPE SUPPORT", funcCall), expected));
    }
}

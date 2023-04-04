package grok_connect.providers.arguments_provider;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.DataFrameBuilder;
import grok_connect.providers.utils.DateParser;
import grok_connect.providers.utils.FuncCallBuilder;
import grok_connect.providers.utils.Parser;
import org.junit.jupiter.api.Named;
import org.junit.jupiter.params.provider.Arguments;
import serialization.BigIntColumn;
import serialization.BoolColumn;
import serialization.DataFrame;
import serialization.DateTimeColumn;
import serialization.FloatColumn;
import serialization.IntColumn;
import serialization.StringColumn;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

public class ClickHouseObjectsMother {
    private static final Parser parser = new DateParser();

    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new StringColumn(new String[] {"INFORMATION_SCHEMA", "default",
                                "information_schema", "system"}),
                        "table_schema")
                .build();
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> getSchema_ok() {
        String firstColumnName = "table_schema";
        String secondColumnName = "table_name";
        String thirdColumnName = "column_name";
        String fourthColumnName = "data_type";
        String fifthColumnName = "is_view";
        String schema = "default";
        String table = "mock_data";
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(10)
                .setColumn(new StringColumn(), firstColumnName, new String[] {schema, schema,
                        schema, schema, schema, schema, schema,
                        schema, schema, schema})
                .setColumn(new StringColumn(), secondColumnName, new String[] {table, table,
                        table, table, table, table,
                        table, table, table, table})
                .setColumn(new StringColumn(), thirdColumnName, new String[] {"id", "first_name", "last_name", "email",
                        "gender", "ip_address", "bool", "country", "date", "some_number"})
                .setColumn(new StringColumn(), fourthColumnName, new String[] {"UInt64", "String",
                        "String", "String", "String", "String",
                        "Bool", "String", "Date", "Float64"})
                .setColumn(new IntColumn(), fifthColumnName, new Integer[] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0})
                .build();
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> checkListParameterSupport_ok() {
        // list<string>
        List<String> values = new ArrayList<>();
        values.add("Poland");
        values.add("Brazil");
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new BigIntColumn(new String[]{"2", "5", "20"}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Nicholle", "Mitchell", "Lucius",}), "first_name")
                .setColumn(new StringColumn(new String[]{"Karoly", "Haglington", "Edelmann"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"nkaroly1@alexa.com", "mhaglington4@indiegogo.com",
                        "ledelmannj@bravesites.com"}), "email")
                .setColumn(new StringColumn(new String[]{"Female", "Male", "Male"}), "gender")
                .setColumn(new StringColumn(new String[]{"255.233.247.118/32", "209.93.181.190/32", "66.174.30.225/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{false, true, false}), "bool")
                .setColumn(new StringColumn(new String[]{"Poland", "Poland", "Brazil"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles("yyyy-MM-dd", "2014-02-27",
                                "2020-10-09","1999-06-22")),
                        "date")
                .setColumn(new FloatColumn(new Float[]{864.09f, 15.22f, 378.73f}), "some_number")
                .build();
        FuncCall funcCall1 = FuncCallBuilder.getBuilder()
                .addQuery("--input: list<string> values\n" +
                        "SELECT * FROM mock_data WHERE country IN (@values)")
                .addFuncParam("list", "string", "values", values, "")
                .addFuncCallOptionsPattern("country", "", "",
                        null, null, "Poland", "Brazil")
                .build();
        FuncCall funcCall2 = FuncCallBuilder.getBuilder()
                .addQuery("--input: list<string> values = ['Poland','Brazil']\n" +
                        "SELECT * FROM mock_data WHERE country IN (@values)")
                .addFuncParam("list", "string", "values", values, "")
                .addFuncCallOptionsPattern("country", "", "",
                        null, null, "Poland", "Brazil")
                .build();
        return Stream.of(
                Arguments.of(Named.of("type: list<string>; operator: none; pattern: none", funcCall1), expected),
                Arguments.of(Named.of("type: list<string>; operator: none; pattern: none", funcCall2), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_uuidType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[] {"417ddc5d-e556-4d27-95dd-a34d84e46a50"}), "uuid")
                .build();
        return Stream.of(Arguments.of(Named.of("UUID TYPE",
                        FuncCallBuilder.fromQuery("SELECT * FROM uuid_type;")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_numericType_ok() {
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new IntColumn(new Integer[]{127, -128}), "int8_type")
                .setColumn(new IntColumn(new Integer[]{32767, -32768}), "int16_type")
                .setColumn(new IntColumn(new Integer[]{2147483647, -2147483648}), "int32_type")
                .setColumn(new BigIntColumn(new String[] {"9223372036854775807", "-9223372036854775808"}), "int64_type")
                .setColumn(new BigIntColumn(new String[] {"170141183460469231731687303715884105727",
                        "-170141183460469231731687303715884105728"}), "int128_type")
                .setColumn(new BigIntColumn(new String[] {"57896044618658097711785492504343953926634992332820282019728792003956564819967",
                        "-57896044618658097711785492504343953926634992332820282019728792003956564819968"}), "int256_type")
                .build();
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new IntColumn(new Integer[]{255}), "uint8_type")
                .setColumn(new IntColumn(new Integer[]{65535}), "uint16_type")
                .setColumn(new BigIntColumn(new String[]{"4294967295"}), "uint32_type")
                .setColumn(new BigIntColumn(new String[] {"18446744073709551615"}), "uint64_type")
                .setColumn(new BigIntColumn(new String[] {"340282366920938463463374607431768211455"}), "uint128_type")
                .setColumn(new BigIntColumn(new String[] {"11579208923731619542357098500868790785"
                        + "3269984665640564039457584007913129639935"}), "uint256_type")
                .build();
        DataFrame expected3 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new FloatColumn(new Float[]{Float.POSITIVE_INFINITY, Float.NEGATIVE_INFINITY}), "float32_type")
                .setColumn(new FloatColumn(new Float[]{Float.POSITIVE_INFINITY, Float.NaN}), "float64_type")
                .setColumn(new FloatColumn(new Float[]{999.9999f, -999.9999f}), "decimal32_type")
                .setColumn(new FloatColumn(new Float[]{1.0E9f, -1.0E9f}), "decimal64_type")
                .setColumn(new FloatColumn(new Float[]{2.11111117E11f, -2.11111117E10f}), "decimal128_type")
                .setColumn(new FloatColumn(new Float[]{9.9999998E17f, -1.0E21f}), "decimal256_type")
                .build();
        return Stream.of(Arguments.of(Named.of("SIGNED INTEGER TYPES",
                FuncCallBuilder.fromQuery("SELECT * FROM SIGNED_INTEGER_TYPES ORDER BY int8_type DESC;")), expected1),
                Arguments.of(Named.of("UNSIGNED INTEGER TYPES",
                        FuncCallBuilder.fromQuery("SELECT * FROM UNSIGNED_INTEGER_TYPES;")), expected2),
                Arguments.of(Named.of("FLOAT TYPES",
                        FuncCallBuilder.fromQuery("SELECT * FROM FLOAT_TYPES ORDER BY float32_type DESC;")), expected3));
    }

    public static Stream<Arguments> checkOutputDataFrame_arrayType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[] {"[1, 2, 3, 4]"}), "array_type")
                .setColumn(new StringColumn(new String[] {"[[24], [421, 12, 4], []]"}), "array_array_type")
                .build();
        return Stream.of(Arguments.of(Named.of("ARRAY TYPE",
                FuncCallBuilder.fromQuery("SELECT * FROM ARRAY_TYPE;")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_tupleType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[] {"[1244, Datagrok]"}), "tuple")
                .build();
        return Stream.of(Arguments.of(Named.of("TUPLE TYPE",
                FuncCallBuilder.fromQuery("SELECT * FROM TUPLE_TYPE;")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_mapType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[] {"{key1=1, key2=10}",
                        "{key1=2, key2=20}", "{key1=3, key2=30}"}), "a")
                .build();
        return Stream.of(Arguments.of(Named.of("MAP TYPE",
                FuncCallBuilder.fromQuery("SELECT * FROM MAP_TYPE;")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_dateTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(new Double[]{
                        parser.parseDateToDouble("yyyy-MM-dd", "2016-06-15")}),
                        "date_type")
                .setColumn(new DateTimeColumn(new Double[]{
                                parser.parseDateToDouble("yyyy-MM-dd", "2016-06-15")}),
                        "date32_type")
                .setColumn(new DateTimeColumn(new Double[]{
                                parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss",
                                        "2016-06-15 23:00:00")}),
                        "datetime_type")
                .setColumn(new DateTimeColumn(new Double[]{
                                parser.parseDateToDouble("yyyy-MM-dd'T'HH:mm:ssX",
                                        "2019-01-01T02:00:00+0200")}),
                        "datetime64_type")
                .build();
        return Stream.of(Arguments.of(Named.of("DATE TYPES",
                FuncCallBuilder.fromQuery("SELECT * FROM DATE_TYPES;")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_nestedType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd",
                                "2016-01-01")}),
                        "EventDate")
                .setColumn(new BigIntColumn(new String[]{"123"}), "UserID")
                .setColumn(new StringColumn(new String[]{"[price, color]"}), "Attrs.Key")
                .setColumn(new StringColumn(new String[]{"[high, red]"}), "Attrs.Value")
                .build();
        return Stream.of(Arguments.of(Named.of("NESTED TYPE",
                FuncCallBuilder.fromQuery("SELECT * FROM NESTED_TYPE;")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_geoTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[]{"[10.0, 10.0]"}), "p")
                .setColumn(new StringColumn(new String[]{"[[0.0, 0.0][10.0, 0.0][10.0, 10.0][0.0, 10.0]]"}), "r")
                .setColumn(new StringColumn(new String[]{"[[[20.0, 20.0][50.0, 20.0][50.0, 50.0][20.0, 50.0]]"
                        + "[[30.0, 30.0][50.0, 50.0][50.0, 30.0]]]"}), "pg")
                .setColumn(new StringColumn(new String[]{"[[[[0.0, 0.0][10.0, 0.0][10.0, 10.0][0.0, 10.0]]]"
                                + "[[[20.0, 20.0][50.0, 20.0][50.0, 50.0][20.0, 50.0]][[30.0, 30.0][50.0, 50.0][50.0, 30.0]]]]"}),
                        "mpg")
                .build();
        return Stream.of(Arguments.of(Named.of("GEO TYPE",
                FuncCallBuilder.fromQuery("SELECT * FROM GEO_TYPES;")), expected));
    }
}

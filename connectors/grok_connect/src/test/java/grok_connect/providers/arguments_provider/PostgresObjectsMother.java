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
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Base64;
import java.util.stream.Stream;

/**
 * Provides data for Postgres provider tests
 */
public class PostgresObjectsMother {
    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[] {"public", "pg_catalog",
                        "information_schema"}), "table_schema")
                .build();
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> getSchema_ok() {
        // based on test script postgres_basic_types.sql
        String firstColumnName = "table_schema";
        String secondColumnName = "table_name";
        String thirdColumnName = "column_name";
        String fourthColumnName = "data_type";
        String fifthColumnName = "is_view";
        String catalog = "datagrok";
        String schema = "public";
        String table = "mock_data";
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(10)
                .setColumn(new StringColumn(), "table_catalog", new String[] {catalog, catalog,
                        catalog, catalog, catalog, catalog, catalog,
                        catalog, catalog, catalog})
                .setColumn(new StringColumn(), firstColumnName, new String[] {schema, schema,
                        schema, schema, schema, schema, schema,
                        schema, schema, schema})
                .setColumn(new StringColumn(), secondColumnName, new String[] {table, table,
                        table, table, table, table,
                        table, table, table, table})
                .setColumn(new StringColumn(), thirdColumnName, new String[] {"id", "first_name", "last_name", "email",
                        "gender", "ip_address", "bool", "country", "date", "some_number"})
                .setColumn(new StringColumn(), fourthColumnName, new String[] {"bigint", "character varying",
                        "character varying", "character varying", "character varying", "cidr",
                        "boolean", "character varying", "date", "numeric"})
                .setColumn(new IntColumn(), fifthColumnName, new Integer[] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0})
                .build();
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_arrayType_ok() {
        // based on test script postgres_array.sql
        String firstColumnName = "name";
        String secondColumnName = "pay_by_quarter";
        String thirdColumnName = "schedule";
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new StringColumn(new String[]{"Bill", "Carol"}), firstColumnName)
                .setColumn(new StringColumn(new String[] {"{10000,10000,10000,10000}",
                        "{20000,25000,25000,25000}"}), secondColumnName)
                .setColumn(new StringColumn(new String[] {"{{meeting,lunch},{training,presentation}}",
                        "{{breakfast,consulting},{meeting,lunch}}"}), thirdColumnName)
                .build();
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new IntColumn(new Integer[] {10000, 25000}), secondColumnName)
                .build();
        DataFrame expected3 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[] {"{{lunch},{presentation}}"}), thirdColumnName)
                .build();
        return Stream.of(
                Arguments.of(Named.of("ARRAY TYPE SUPPORT",
                        FuncCallBuilder.fromQuery("SELECT * FROM sal_emp;")), expected1),
                Arguments.of(Named.of("ARRAY'S INTERNAL TYPE SUPPORT#1",
                        FuncCallBuilder.fromQuery("SELECT pay_by_quarter[3] FROM sal_emp;")), expected2),
                Arguments.of(Named.of("ARRAY'S INTERNAL TYPE SUPPORT#2",
                        FuncCallBuilder.fromQuery("SELECT schedule[:2][2:] FROM sal_emp WHERE "
                        + "name = 'Bill';")), expected3)
        );
    }

    public static Stream<Arguments> checkOutputDataFrame_basicTypes_ok() {
        Parser parser = new DateParser();
        String datePattern = "yyyy-MM-dd";
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(10)
                .setColumn(new BigIntColumn(
                                Stream.iterate(1, x -> x + 1).limit(10)
                                        .map(String::valueOf).toArray(String[]::new)),
                        "id")
                .setColumn(new StringColumn(new String[]{"Burk", "Nicholle", "Orlando", "Gothart",
                        "Mitchell", "Jeromy", "Joela", "Darren", "Marlie", "Scottie"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Kemery", "Karoly", "Westgate",
                                "Cokayne", "Haglington", "Twinn", "Cornau", "Juares", "Mayze", "Formilli"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"bkemery0@businesswire.com", "nkaroly1@alexa.com",
                        "owestgate2@dedecms.com", "gcokayne3@plala.or.jp", "mhaglington4@indiegogo.com",
                        "jtwinn5@globo.com", "jcornau6@imgur.com", "djuares7@hexun.com",
                        "mmayze8@google.com.au", "sformilli9@aol.com"}), "email")
                .setColumn(new StringColumn(new String[]{"Male", "Female", "Polygender", "Male", "Male", "Male",
                        "Female", "Male", "Female", "Male"}), "gender")
                .setColumn(new StringColumn(new String[]{"249.64.22.121/32", "255.233.247.118/32",
                                "75.0.252.254/32", "196.83.12.163/32", "209.93.181.190/32", "25.13.2.132/32",
                                "195.47.88.236/32", "94.170.16.96/32", "68.41.25.65/32", "101.241.191.228/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{true, false, false, true, true, true,
                        false, false, false, false}), "bool")
                .setColumn(new StringColumn(new String[]{"China", "Poland", "Netherlands",
                        "Philippines", "Poland", "Serbia", "Indonesia", "China",
                        "France", "Vietnam"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern,"2017-09-20",
                        "2014-02-27", "2020-09-03", "2001-01-31", "2020-10-09",
                        "2014-10-04", "2020-03-19", "2011-04-09", "2011-11-10", "2003-01-04")), "date")
                .setColumn(new FloatColumn(new Float[]{510.32f, 864.09f, 822.7f, 251.05f, 15.22f,
                        378.4f, 349.11f, 631.89f, 561.72f, 978.01f}), "some_number")
                .build();
        return Stream.of(Arguments.of(Named.of("VARCHAR, BOOLEAN, BIGINT, DATE, CIDR, NUMERIC TYPES SUPPORT",
                FuncCallBuilder.fromQuery("SELECT * FROM mock_data LIMIT 10;")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_bitStringType_ok() {
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new BigIntColumn(new String[] {"101", "001"}), "a")
                .setColumn(new BigIntColumn(new String[] {"0011", "101"}), "b")
                .build();
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new BoolColumn(new Boolean[] {true, false}), "c")
                .setColumn(new BigIntColumn(new String[] {"101", "0"}), "d")
                .build();
        return Stream.of(Arguments.of(Named.of("BIT STRING WITH FIXED(>1) AND VARYING LENGTH SUPPORT",
                        FuncCallBuilder.fromQuery("SELECT * FROM test1;")), expected1),
                Arguments.of(Named.of("BIT STRING WITH FIXED(1) AND VARYING LENGTH SUPPORT",
                        FuncCallBuilder.fromQuery("SELECT * FROM test2")), expected2));
    }

    public static Stream<Arguments> checkOutputDataFrame_byteAType_ok() {
        Path path = Paths.get("src/test/resources/scripts/postgres/file.txt");
        byte[] file;
        try {
            file = Files.readAllBytes(path);
        } catch (IOException e) {
            throw new RuntimeException("Something went wrong when reading file " + path, e);
        }
        String data = Base64.getEncoder().encodeToString(file);// or if we know that it is text file,
        // we can use new String(data, StandardCharset.UTF-8)
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[] {data}), "data")
                .build();
        return Stream.of(Arguments.of(Named.of("BYTEA TYPE SUPPORT",
                FuncCallBuilder.fromQuery("SELECT encode(data, 'base64') as data FROM bytea_data;")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_compositeType_ok() {
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[]{"(\"fuzzy dice\",42,1.99)"}), "item")
                .setColumn(new IntColumn(new Integer[]{1000}), "count")
                .build();
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[]{"fuzzy dice"}), "name")
                .setColumn(new IntColumn(new Integer[]{42}), "supplier_id")
                .setColumn(new FloatColumn(new Float[]{1.99f}), "price")
                .setColumn(new IntColumn(new Integer[]{1000}), "count")
                .build();
        return Stream.of(
                Arguments.of(Named.of("COMPOSITE TYPE SUPPORT",
                        FuncCallBuilder.fromQuery("SELECT * FROM on_hand;")), expected1),
                Arguments.of(Named.of("COMPOSITE TYPE FIELDS TYPES SUPPORT",
                        FuncCallBuilder.fromQuery("SELECT (item).name, (item).supplier_id, "
                        + "(item).price, count FROM on_hand")), expected2)
        );
    }

    public static Stream<Arguments> checkOutputDataFrame_dateTypes_ok() {
        Parser parser = new DateParser();
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd",
                        "1999-01-08")}), "date")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("hh:mm:ss.SSS",
                        "04:05:06.789")}), "time")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss",
                        "1999-01-08 04:05:06"),}), "stamp")
                .setColumn(new StringColumn(new String[] {"1 years 5 mons 5 days 0 hours 0 mins 0.0 secs"}), "interval")
                .build();
        return Stream.of(Arguments.of(Named.of("DATE, TIME, TIMESTAMP, INTERVAL TYPES SUPPORT",
                FuncCallBuilder.fromQuery("SELECT * FROM dates")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_jsonbType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(
                        new String[]{"{\"phones\": [{\"type\": \"mobile\", \"phone\": \"001001\"}, {\"type\": \"fix\", \"phone\": \"002002\"}]}",
                                "{\"bar\": \"baz\", \"active\": false, \"balance\": 7.77}", "{\"reading\": 0.00001230}"}), "data")
                .build();
        return Stream.of(Arguments.of(Named.of("JSONB TYPE SUPPORT",
                FuncCallBuilder.fromQuery("SELECT * FROM jsonb_data")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_numericType_ok() {
        String s1 = "-9223372036854775808";
        String s2 = "9223372036854775807";
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new FloatColumn(new Float[]{Float.NaN, 2.22222222222f, 222222222222222222.2f}),
                        "data")
                .build();
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(4)
                .setColumn(new FloatColumn(new Float[]{Float.NEGATIVE_INFINITY, Float.POSITIVE_INFINITY,
                        0.0f, Float.POSITIVE_INFINITY}), "data")
                .build();
        DataFrame expected3 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new FloatColumn(new Float[]{1E-37f, 1E+37f}), "data")
                .build();
        DataFrame expected4 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new BigIntColumn(new String[]{s1,
                                s2}),
                        "data")
                .build();
        DataFrame expected5 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new FloatColumn(new Float[] {2f, 2222f}), "data")
                .build();
        DataFrame expected6 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new FloatColumn(new Float[] {22.22f, 55.55f}), "data")
                .build();
        return Stream.of(
                Arguments.of(Named.of("NUMERIC TYPE SUPPORT WITHOUT PRECISION AND SCALE",
                        FuncCallBuilder.fromQuery("SELECT * FROM numeric_data")), expected1),
                Arguments.of(Named.of("DOUBLE TYPE SUPPORT",
                        FuncCallBuilder.fromQuery("SELECT * FROM doubles")), expected2),
                Arguments.of(Named.of("REAL TYPE SUPPORT",
                        FuncCallBuilder.fromQuery("SELECT * FROM reals")), expected3),
                Arguments.of(Named.of("BIGINT TYPE SUPPORT",
                        FuncCallBuilder.fromQuery("SELECT * FROM bigint_data")), expected4),
                Arguments.of(Named.of("NUMERIC TYPE SUPPORT WITH PRECISION",
                        FuncCallBuilder.fromQuery("SELECT * FROM numeric_data_precision")), expected5),
                Arguments.of(Named.of("NUMERIC TYPE SUPPORT WITH PRECISION AND SCALE",
                        FuncCallBuilder.fromQuery("SELECT * FROM numeric_data_precision_scale")), expected6)
        );
    }

    public static Stream<Arguments> checkOutputDataFrame_serialType_ok() {
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new IntColumn(new Integer[]{1, 2}), "id")
                .setColumn(new StringColumn(new String[]{"Fiat", "Honda"}), "brand")
                .build();
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new BigIntColumn(new String[]{"1", "2"}), "id")
                .setColumn(new StringColumn(new String[]{"Fiat", "Honda"}), "brand")
                .build();
        return Stream.of(
                Arguments.of(Named.of("SMALL SERIAL TYPE SUPPORT", FuncCallBuilder.fromQuery("SELECT * FROM cars_small")), expected1),
                Arguments.of(Named.of("SERIAL TYPE SUPPORT", FuncCallBuilder.fromQuery("SELECT * FROM cars")), expected1),
                Arguments.of(Named.of("BIG SERIAL TYPE SUPPORT", FuncCallBuilder.fromQuery("SELECT * FROM cars_big")), expected2)
        );
    }

    public static Stream<Arguments> checkOutputDataFrame_uuidType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[]{"6ecd8c99-4036-403d-bf84-cf8400f67836",
                        "c81d4e2e-bcf2-11e6-869b-7df92533d2db", "237e9877-e79b-12d4-a765-321741963000"}), "id")
                .build();
        return Stream.of(Arguments.of(Named.of("UUID TYPE SUPPORT", FuncCallBuilder.fromQuery("SELECT * FROM uuid_data")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_xmlType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[]{"<foo>Hello World!</foo>",
                        "abc<foo>bar</foo><bar>foo</bar>",
                        "<book><title>Manual</title><chapter>...</chapter></book>"}), "data")
                .build();
        return Stream.of(Arguments.of(
                Named.of("XML TYPE SUPPORT", FuncCallBuilder.fromQuery("SELECT * FROM xml_data")), expected
            )
        );
    }

    public static Stream<Arguments> checkPostgresOperatorsSupport_ok() {
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new FloatColumn(new Float[]{5.0f}), "value")
                .build();
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new BigIntColumn(new String[]{"1"}), "id")
                .setColumn(new StringColumn(new String[]{"{\"reading\": 0.00001230}"}), "json_data")
                .setColumn(new FloatColumn(new Float[]{2.0f}), "length")
                .setColumn(new StringColumn(new String[]{"(0.0,0.0)"}), "center")
                .build();
        FuncCall funcCall2 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string id = \"=1\" {pattern: int}\n"
                        + "SELECT id, json_data, @-@ path_data AS length, @@ circle_data AS center FROM operators WHERE @id(id)\n"
                        + "--end")
                .addFuncParam("string", "", "id", "=1", "int")
                .addFuncCallOptionsPattern("id", "=1", "=",
                        null, null, 1)
                .build();
        DataFrame expected3 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new BigIntColumn(new String[]{"2"}), "id")
                .setColumn(new StringColumn(new String[]{"{\"bar\": \"baz\", \"active\": false, \"balance\": 7.77}"}), "json_data")
                .setColumn(new FloatColumn(new Float[]{2.4122634f}), "length")
                .setColumn(new StringColumn(new String[]{"(0.0,0.0)"}), "center")
                .build();
        FuncCall funcCall3 = FuncCallBuilder.getBuilder()
                .addQuery("--input: int id = 2\n"
                        + "SELECT id, json_data, @-@ path_data AS length, @@ circle_data AS center FROM operators "
                        + "WHERE id = @id AND json_data @> '{\"bar\": \"baz\"}'::jsonb\n"
                        + "--end")
                .addFuncParam("int", "", "id", 2, "")
                .build();
        return Stream.of(
                Arguments.of(Named.of("ABSOLUTE VALUE @", FuncCallBuilder.fromQuery("SELECT @ -5.0 AS VALUE")), expected1),
                Arguments.of(Named.of("PATH LENGTH @-@, CIRCLE CENTER @@, PATTERN @ID(ID)", funcCall2), expected2),
                Arguments.of(Named.of("PATH LENGTH @-@, CIRCLE CENTER @@, JSON @>, ID = @ID", funcCall3), expected3)
        );
    }
}

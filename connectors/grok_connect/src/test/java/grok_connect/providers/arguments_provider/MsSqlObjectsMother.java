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
import java.time.Year;
import java.util.stream.Stream;

public class MsSqlObjectsMother {
    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(13)
                .setColumn(new StringColumn(new String[] {
                                "db_accessadmin",
                                "db_backupoperator",
                                "db_datareader",
                                "db_datawriter",
                                "db_ddladmin",
                                "db_denydatareader",
                                "db_denydatawriter",
                                "db_owner",
                                "db_securityadmin",
                                "dbo",
                                "guest",
                                "INFORMATION_SCHEMA",
                                "sys"}),
                        "SCHEMA_NAME")
                .build();
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> getSchema_ok() {
        String firstColumnName = "table_schema";
        String secondColumnName = "table_name";
        String thirdColumnName = "column_name";
        String fourthColumnName = "data_type";
        String fifthColumnName = "is_view";
        String schema = "dbo";
        String table = "MOCK_DATA";
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(9)
                .setColumn(new StringColumn(), firstColumnName, new String[] {schema, schema,
                        schema, schema, schema, schema,
                        schema, schema, schema})
                .setColumn(new StringColumn(), secondColumnName, new String[] {table, table,
                        table, table, table, table,
                        table, table, table})
                .setColumn(new StringColumn(), thirdColumnName, new String[] {"id", "first_name", "last_name", "email",
                        "gender", "ip_address", "country", "date", "some_number"})
                .setColumn(new StringColumn(), fourthColumnName, new String[] {"bigint", "varchar",
                        "varchar", "varchar", "varchar", "varchar", "varchar", "date", "numeric"})
                .setColumn(new IntColumn(), fifthColumnName, new Integer[] {0, 0, 0, 0, 0, 0, 0, 0, 0})
                .build();
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> checkMultipleParametersSupport_ok() {
        Parser parser = new DateParser();
        String datePattern = "yyyy-MM-dd";
        // --input: string first_name = "starts with p" {pattern: string}
        //--input: string id = ">1" {pattern :int}
        //--input: bool bool = false
        //--input: string email = "contains com" {pattern: string}
        //--input: string some_number = ">20" {pattern: double}
        //--input: string country = "in (Indonesia)" {pattern: string}
        //--input: string date = "before 1/1/2022" {pattern: datetime}
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new BigIntColumn(new String[]{"13"}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Pail"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Boxell"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"pboxellc@moonfruit.com"}), "email")
                .setColumn(new StringColumn(new String[]{"Genderqueer"}), "gender")
                .setColumn(new StringColumn(new String[]{"2.37.160.155/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{false}), "bool")
                .setColumn(new StringColumn(new String[]{"Indonesia"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2012-01-14")),
                        "date")
                .setColumn(new FloatColumn(new Float[]{73.47f}), "some_number")
                .build();
        FuncCall funcCall1 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string first_name = \"starts with p\" {pattern: string}\n"
                        + "--input: string id = \">1\" {pattern :int}\n"
                        + "--input: string email = \"contains com\" {pattern: string}\n"
                        + "--input: string some_number = \">20\" {pattern: double}\n"
                        + "--input: string country = \"in (Indonesia)\" {pattern: string}\n"
                        + "--input: string date = \"before 1/1/2022\" {pattern: datetime}\n"
                        + "SELECT * FROM mock_data WHERE @first_name(first_name) AND @id(id) "
                        + "AND @email(email) AND @some_number(some_number) "
                        + "AND @country(country) AND @date(date)\n"
                        + "--end")
                .addFuncParam("string", "", "first_name", "starts with p", "string")
                .addFuncParam("string", "","id", ">1", "string")
                .addFuncParam("string", "","email", "contains com", "string")
                .addFuncParam("string", "","some_number", ">20", "double")
                .addFuncParam("string", "","country", "in (Indonesia)", "string")
                .addFuncParam("string", "","date", "before 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("first_name", "starts with p",
                        "starts with", null, null, "p")
                .addFuncCallOptionsPattern("id", ">1", ">", null,
                        null, "1")
                .addFuncCallOptionsPattern("email", "contains com",
                        "contains", null, null, "com")
                .addFuncCallOptionsPattern("some_number", ">20", ">", null,
                        null, 20)
                .addFuncCallOptionsPattern("country", "in (Indonesia)", "in",
                        null, null, "Indonesia")
                .addFuncCallOptionsPattern("date", "before 1/1/2022", "before",
                        true, true, Year.of(2022).atMonth(1).atDay(1).toString())
                .build();
        return Stream.of(Arguments.of(Named.of("type: multiple; operator: multiple; pattern: multiple", funcCall1),
                expected1));
    }

    public static Stream<Arguments> checkOutputDataFrame_dateTypes_ok() {
        Parser parser = new DateParser();
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd",
                        "1900-01-01")}), "date_data")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss",
                        "1900-01-01 00:00:00")}), "datetime_data")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd'T'HH:mm:ss.SSSSSSS",
                        "2007-05-02T19:58:47.123"),}), "datetime_data2")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("HH:mm:ss.SS",
                        "00:00:00.00")}), "time_data")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss X",
                        "1900-01-01 00:00:00 +02:00"),}), "datetimeoffset_data")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss",
                        "1955-12-13 12:43:00")}), "smalldatetime_data")
                .build();
        return Stream.of(Arguments.of(Named.of("DATE, DATETIME, DATETIME2, TIME, DATETIMEOFFSET, SMALL DATETIME " +
                        "TYPES SUPPORT",
                FuncCallBuilder.fromQuery("SELECT * FROM date_types")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_characterTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[]{"Datagrok            "}), "char_data")
                .setColumn(new StringColumn(new String[]{"Hello, World!"}), "varchar_data")
                .setColumn(new StringColumn(new String[]{"Lorem ipsum"}), "text_data")
                .setColumn(new StringColumn(new String[]{"Hello, World!   "}), "nchar_data")
                .setColumn(new StringColumn(new String[]{"Hello, Datagrok!"}), "nvarchar_data")
                .setColumn(new StringColumn(new String[]{"Hello, Datagrok!"}), "ntext_data")
                .build();
        return Stream.of(Arguments.of(Named.of("CHARACTER TYPES SUPPORT",
                FuncCallBuilder.fromQuery("SELECT * FROM character_types")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_xmlType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new StringColumn(new String[]{"<foo>Hello World!</foo>",
                        "<book><title>Manual</title><chapter>...</chapter></book>"}),
                        "xml_data")
                .build();
        return Stream.of(Arguments.of(Named.of("XML TYPE SUPPORT",
                FuncCallBuilder.fromQuery("SELECT * FROM xml_type")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_binaryTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[]{"Datagrok"}),
                        "binary_data")
                .setColumn(new IntColumn(new Integer[]{123456}),
                        "varbinary_data")
                .build();
        return Stream.of(Arguments.of(Named.of("BINARY TYPES SUPPORT",
                FuncCallBuilder.fromQuery("SELECT CONVERT(VARCHAR(max), binary_data, 0) as binary_data, "
                        + "CONVERT(int, varbinary_data, 1) as varbinary_data FROM binary_types")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_numericTypes_ok() {
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new BigIntColumn(new String[]{"9223372036854775807"}), "bigint_data")
                .setColumn(new IntColumn(new Integer[]{2147483647}), "int_data")
                .setColumn(new IntColumn(new Integer[]{32767}), "smallint_data")
                .setColumn(new IntColumn(new Integer[]{123}), "tinyint_data")
                .setColumn(new BoolColumn(new Boolean[]{true}), "bit_data")
                .setColumn(new FloatColumn(new Float[]{123.22f}), "decimal_data")
                .setColumn(new FloatColumn(new Float[]{12345.12000f}), "numeric_data")
                .build();
        FuncCall funcCall1 = FuncCallBuilder.fromQuery("SELECT * FROM numeric_types");
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new FloatColumn(new Float[]{-1.79E30f, 34636.34661f }), "float_data1")
                .setColumn(new FloatColumn(new Float[]{-0.0f, 0.0f}), "float_data2")
                .setColumn(new FloatColumn(new Float[]{124124.23555f, 0.0f}), "real_data")
                .build();
        FuncCall funcCall2 = FuncCallBuilder.fromQuery("SELECT * FROM float_types");
        return Stream.of(Arguments.of(Named.of("NUMERIC TYPES SUPPORT", funcCall1), expected1),
                Arguments.of(Named.of("FLOAT TYPES SUPPORT", funcCall2), expected2));
    }

    public static Stream<Arguments> checkOutputDataFrame_moneyTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new FloatColumn(new Float[]{922337203685477.5807f, -922337203685477.5808f, 346.46f}), "money_type")
                .setColumn(new FloatColumn(new Float[]{214748.3647f, -214748.3648f, 160.0f}), "small_money")
                .build();
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM money_types");
        return Stream.of(Arguments.of(Named.of("MONEY TYPES SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_spatialTypes_ok() {
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new IntColumn(new Integer[]{1, 2}), "id")
                .setColumn(new StringColumn(new String[]{"LINESTRING (100 100, 20 180, 180 180)",
                        "POLYGON ((0 0, 150 0, 150 150, 0 150, 0 0))"}), "GeomCol1")
                .setColumn(new StringColumn(new String[]{"LINESTRING (100 100, 20 180, 180 180)",
                        "POLYGON ((0 0, 150 0, 150 150, 0 150, 0 0))"}), "GeomCol2")
                .build();
        FuncCall funcCall1 = FuncCallBuilder.fromQuery("SELECT id, GeomCol1.ToString() as GeomCol1, "
                + "GeomCol2 FROM SpatialTable1");
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new IntColumn(new Integer[]{1, 2}), "id")
                .setColumn(new StringColumn(new String[]{"LINESTRING (-122.36 47.656, -122.343 47.656)",
                        "POLYGON ((-122.358 47.653, -122.348 47.649, -122.348 47.658, -122.358 47.658, -122.358 47.653))"}),
                        "GeogCol1")
                .setColumn(new StringColumn(new String[]{"LINESTRING (-122.36 47.656, -122.343 47.656)",
                                "POLYGON ((-122.358 47.653, -122.348 47.649, -122.348 47.658, -122.358 47.658, -122.358 47.653))"}),
                        "GeogCol2")
                .build();
        FuncCall funcCall2 = FuncCallBuilder.fromQuery("SELECT id, GeogCol1.ToString() as GeogCol1, "
                + "GeogCol2 FROM SpatialTable2");
        return Stream.of(Arguments.of(Named.of("GEOMETRY TYPE SUPPORT", funcCall1), expected1),
                Arguments.of(Named.of("GEOGRAPHY TYPE SUPPORT", funcCall2), expected2));
    }
}

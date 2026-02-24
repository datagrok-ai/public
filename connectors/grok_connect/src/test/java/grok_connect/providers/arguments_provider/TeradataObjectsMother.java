package grok_connect.providers.arguments_provider;

import grok_connect.connectors_info.FuncCall;
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
import java.time.LocalDate;
import java.time.Year;
import java.time.temporal.TemporalAdjusters;
import java.util.stream.Stream;

@SuppressWarnings("unused")
public class TeradataObjectsMother {
    private static final Parser parser = new DateParser();

    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn("table_schema", new String[] {"tdwm", "TDStats", "tapidb",
                                "TDMaps", "val", "TDBCMgmt", "SYSLIB", "mldb", "TD_ANALYTICS_DB",
                                "TEST", "TD_SYSXML", "SystemFe", "dbcmngr", "TD_SERVER_DB", "SYSUDTLIB",
                                "TDQCD", "LockLogShredder", "SYSUIF", "Sys_Calendar", "TD_SYSFNLIB",
                                "TDaaS_DB", "DBC", "SYSSPATIAL", "SQLJ", "SYSJDBC", "SYSBAR", "SysAdmin"}));
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> getSchema_ok() {
        String firstColumnName = "table_schema";
        String secondColumnName = "table_name";
        String thirdColumnName = "column_name";
        String fourthColumnName = "data_type";
        String fifthColumnName = "is_view";
        String schema = "TEST";
        String table = "MOCK_DATA";
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn(firstColumnName, new String[] {schema, schema,
                        schema, schema, schema, schema, schema,
                        schema, schema}),
                new StringColumn(secondColumnName, new String[] {table, table,
                        table, table, table, table,
                        table, table, table}),
                new StringColumn(thirdColumnName, new String[] {"first_name", "ip_address",
                        "last_name", "some_number", "gender", "email", "country", "id", "dat"}),
                new StringColumn(fourthColumnName, new String[] {"VARCHAR(50)", "VARCHAR(50)",
                        "VARCHAR(50)", "NUMBER(5,2)", "VARCHAR(50)", "VARCHAR(50)",
                        "VARCHAR(50)", "BIGINT", "DATE"}),
                new IntColumn(fifthColumnName, new Integer[] {0, 0, 0, 0, 0, 0, 0, 0, 0}));
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> checkDatesParameterSupport_ok() {
        String datePattern = "yyyy-MM-dd";
        LocalDate now = LocalDate.now();
        int dayOfWeek = now.getDayOfWeek().getValue();
        int dayOfMonth = now.getDayOfMonth();
        int dayOfYear = now.getDayOfYear();
        LocalDate firstDayOfWeek = now.minusDays(dayOfWeek - 1);
        LocalDate lastDayOfWeek = now.plusDays(7 - dayOfWeek);
        LocalDate firstDayOfMonth = now.with(TemporalAdjusters.firstDayOfMonth());
        LocalDate lastDayOfMonth = now.with(TemporalAdjusters.lastDayOfMonth());
        LocalDate firstDayOfYear = now.with(TemporalAdjusters.firstDayOfYear());
        LocalDate lastDayOfYear = now.with(TemporalAdjusters.lastDayOfYear());
        LocalDate yesterday = now.minusDays(1);
        LocalDate dayOfLastYear = now.minusDays(150);
        // --input: string date = "today" {pattern: datetime}
        DataFrame expected1 = DataFrame.fromColumns(
                new DateTimeColumn("dat", new Double[]{parser.parseDateToDouble(datePattern, now.toString())}));
        FuncCall funcCall1 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"today\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "", "dat", "today", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        now.toString(), now.plusDays(1).toString())
                .build();
        // --input: string date = "this week" {pattern: datetime}
        DataFrame expected2 = DataFrame.fromColumns(
                new DateTimeColumn("dat", parser.parseDatesToDoubles(datePattern,
                                dayOfWeek == 1 ? null : yesterday.toString(),
                                now.toString(),
                                lastDayOfWeek.equals(now) ? null : lastDayOfWeek.toString())));
        FuncCall funcCall2 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"this week\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY dat\n"
                        + "--end")
                .addFuncParam("string","", "dat", "this week", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        firstDayOfWeek.toString(),
                        lastDayOfWeek.plusDays(1).toString())
                .build();
        // --input: string date = "this month" {pattern: datetime}
        DataFrame expected3 = DataFrame.fromColumns(
                new DateTimeColumn("dat", parser.parseDatesToDoubles(datePattern,
                                dayOfMonth == 1 ? null : yesterday.toString(),
                                now.toString(),
                                lastDayOfWeek.getMonthValue() >  lastDayOfMonth.getMonthValue() || lastDayOfWeek.equals(now)?
                                        null : lastDayOfWeek.toString())));
        FuncCall funcCall3 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"this month\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY dat\n"
                        + "--end")
                .addFuncParam("string", "","dat", "this month", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        firstDayOfMonth.toString(),
                        lastDayOfMonth.plusDays(1).toString())
                .build();
        // --input: string date = "this year" {pattern: datetime}
        DataFrame expected4 = DataFrame.fromColumns(
                new DateTimeColumn("dat", parser.parseDatesToDoubles(datePattern,
                                dayOfLastYear.getYear() == now.getYear() ? dayOfLastYear.toString() : null,
                                dayOfYear == 1 ? null : yesterday.toString(),
                                now.toString(),
                                lastDayOfWeek.getYear() >  now.getYear() || lastDayOfWeek.equals(now)?
                                        null : lastDayOfWeek.toString())));
        FuncCall funcCall4 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"this year\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY dat\n"
                        + "--end")
                .addFuncParam("string", "","dat", "this year", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        firstDayOfYear.toString(),
                        lastDayOfYear.plusDays(1).toString())
                .build();
        // --input: string date = "yesterday" {pattern: datetime}
        DataFrame expected5 = DataFrame.fromColumns(
                new DateTimeColumn("dat", parser.parseDatesToDoubles(datePattern, yesterday.toString())));
        FuncCall funcCall5 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string date = \"yesterday\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY dat\n"
                        + "--end")
                .addFuncParam("string", "","dat", "yesterday", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        yesterday.toString(),
                        now.toString())
                .build();
        // --input: string date = "last year" {pattern: datetime}
        DataFrame expected6 = DataFrame.fromColumns(
                new DateTimeColumn("dat", parser.parseDatesToDoubles(datePattern,
                                dayOfLastYear.getYear() < now.getYear() ? dayOfLastYear.toString() : null,
                                yesterday.getYear() < now.getYear() ? yesterday.toString() : null)));
        FuncCall funcCall6 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"last year\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY dat\n"
                        + "--end")
                .addFuncParam("string", "", "dat", "last year", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        firstDayOfYear.minusYears(1).toString(), firstDayOfYear.toString())
                .build();
        // --input: string date = "anytime" {pattern: datetime}
        DataFrame expected7 = DataFrame.fromColumns(
                new DateTimeColumn("dat", parser.parseDatesToDoubles(datePattern, "2021-04-09",
                                dayOfLastYear.toString(),
                                yesterday.toString(),
                                now.toString(),
                                lastDayOfWeek.equals(now) ? null : lastDayOfWeek.toString()
                                )));
        FuncCall funcCall7 = FuncCallBuilder.getBuilder()
                .addQuery( "--input: string dat = \"anytime\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY dat\n"
                        + "--end")
                .addFuncParam("string", "", "dat", "anytime", "datetime")
                .addFuncCallOptionsPattern("dat", "", "none", true, true)
                .build();
        // --input: string date = "2021-2022" {pattern: datetime}
        DataFrame expected8 = DataFrame.fromColumns(
                new DateTimeColumn("dat", parser.parseDatesToDoubles(datePattern, "2021-04-09")));

        FuncCall funcCall8 = FuncCallBuilder.getBuilder()
                .addQuery( "--input: string dat = \"2021-2021\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY dat\n"
                        + "--end")
                .addFuncParam("string", "", "dat", "2021-2022", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, true,
                        Year.of(2021).atDay(1).toString(),
                        Year.of(2022).atDay(1).toString())
                .build();
        // --input: string date = "before 1/1/2022" {pattern: datetime}

        FuncCall funcCall9 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"before 1/1/2022\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY dat\n"
                        + "--end")
                .addFuncParam("string", "","dat", "before 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("dat", "", "before", true, true,
                        Year.of(2022).atDay(1).toString())
                .build();
        // --input: string date = "after 1/1/2022" {pattern: datetime}
        DataFrame expected9 = DataFrame.fromColumns(
                new DateTimeColumn("dat", parser.parseDatesToDoubles(datePattern,
                                dayOfLastYear.toString(),
                                yesterday.toString(),
                                now.toString(),
                                lastDayOfWeek.equals(now) ? null : lastDayOfWeek.toString())));
        FuncCall funcCall10 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"after 1/1/2022\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY dat\n"
                        + "--end")
                .addFuncParam("string", "","dat", "after 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("dat", "", "after", true, true,
                        LocalDate.parse("2022-01-01").toString())
                .build();
        // --input: string date = "April 2021" {pattern: datetime}
        FuncCall funcCall11 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"April 2021\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY dat\n"
                        + "--end")
                .addFuncParam("string", "","dat", "April 2021", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        Year.of(2021).atMonth(4).atDay(1).toString(),
                        Year.of(2021).atMonth(5).atDay(1).toString())
                .build();
        return Stream.of(
                Arguments.of(Named.of("type: string; operator: today; pattern: datetime", funcCall1), expected1),
                Arguments.of(Named.of("type: string; operator: this week; pattern: datetime", funcCall2), expected2),
                Arguments.of(Named.of("type: string; operator: this month; pattern: datetime", funcCall3), expected3),
                Arguments.of(Named.of("type: string; operator: this year; pattern: datetime", funcCall4), expected4),
                Arguments.of(Named.of("type: string; operator: yesterday; pattern: datetime", funcCall5), expected5),
                Arguments.of(Named.of("type: string; operator: last year; pattern: datetime", funcCall6), expected6),
                Arguments.of(Named.of("type: string; operator: anytime; pattern: datetime", funcCall7), expected7),
                Arguments.of(Named.of("type: string; operator: range -; pattern: datetime", funcCall8), expected8),
                Arguments.of(Named.of("type: string; operator: before; pattern: datetime", funcCall9), expected8),
                Arguments.of(Named.of("type: string; operator: after; pattern: datetime", funcCall10), expected9),
                Arguments.of(Named.of("type: string; operator: April 2021; pattern: datetime", funcCall11), expected8)
        );
    }

    public static Stream<Arguments> checkArrayTypeSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"1", "2"}),
                new StringColumn("data", new String[]{"ResultArray:INTEGER[121, 255, 1241244]",
                        "ResultArray:INTEGER[0, -124412, 5555]"}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM ARRAY_TYPE");
        return Stream.of(Arguments.of(Named.of("ARRAY TYPE SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkCharacterTypesSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn("ch", new String[]{"Hello     "}),
                new StringColumn("varch", new String[]{"World"}),
                new StringColumn("clb", new String[]{"Datagrok"}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM character_type");
        return Stream.of(Arguments.of(Named.of("CHARACTER TYPES SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkJsonTypeSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"1", "2", "3"}),
                new StringColumn("data", new String[]{"{\"phones\":[{\"type\": \"mobile\", \"phone\": \"001001\"}, "
                        + "{\"type\": \"fix\", \"phone\": \"002002\"}]}",
                        "{\"bar\": \"baz\", \"balance\": 7.77, \"active\":false}", "{\"reading\": 1.230e-5}"}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM JSON_DATA ORDER BY id");
        return Stream.of(Arguments.of(Named.of("JSON TYPE SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkDateTypesSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new DateTimeColumn("dat", new Double[]{parser.parseDateToDouble("yyyy-MM-dd",
                        "2023-01-01")}),
                new DateTimeColumn("time_data", new Double[]{parser.parseDateToDouble("HH:mm:ss.SSS",
                        "12:55:33.333")}),
                new DateTimeColumn("stamp", new Double[]{parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss.SSS",
                        "2023-01-01 12:55:33.333")}),
                new DateTimeColumn("time_zoned", new Double[]{parser.parseDateToDouble("HH:mm:ss.SSSX",
                        "12:55:33.333+02:00")}),
                new DateTimeColumn("zoned_stamp", new Double[]{parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss.SSS",
                        "2000-01-01 11:37:58.222")}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM dates_type");
        return Stream.of(Arguments.of(Named.of("DATE TYPES SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkSpatialTypesSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"1"}),
                new StringColumn("point", new String[]{"POINT (10 20)"}),
                new StringColumn("polygon", new String[]{"POLYGON ((1 1,1 3,6 3,6 0,1 1))"}),
                new StringColumn("linestring", new String[]{"LINESTRING (1 1,2 2,3 3,4 4)"}),
                new StringColumn("multipolygon", new String[]{"MULTIPOLYGON (((1 1,1 3,6 3,6 0,1 1)),((10 5,10 10,20 10,20 5,10 5)))"}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM SPATIAL_TYPES");
        return Stream.of(Arguments.of(Named.of("SPATIAL TYPES SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkXmlTypeSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"1", "2"}),
                new StringColumn("data", new String[]{"<foo>Hello World!</foo>",
                        "<book><title>Manual</title><chapter>...</chapter></book>"}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM XML_DATA");
        return Stream.of(Arguments.of(Named.of("XML TYPE SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkIntegerTypesSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new BigIntColumn("bigint_type", new String[]{"9223372036854775807", "0"}),
                new IntColumn("byte_int", new Integer[] {-128, 127}),
                new IntColumn("small_int", new Integer[] {32760, -23444}),
                new IntColumn("int_type", new Integer[] {2147483647, -2147483647}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM INTEGER_TYPES");
        return Stream.of(Arguments.of(Named.of("INTEGER TYPES SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkFloatTypesSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new FloatColumn("float_type", new Float[]{3.4020234234f, 0.9999f}),
                new FloatColumn("real_type", new Float[]{-3.4028423f, 1.9999f}),
                new FloatColumn("double_type", new Float[]{-0.0043512f, 213515.435545f}),
                new FloatColumn("decimal_type", new Float[]{235234.4646f, 99999999.9999f}),
                new FloatColumn("number_type", new Float[]{999.99f, 0.99f}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM FLOAT_TYPE");
        return Stream.of(Arguments.of(Named.of("FLOAT TYPES SUPPORT", funcCall), expected));
    }
}

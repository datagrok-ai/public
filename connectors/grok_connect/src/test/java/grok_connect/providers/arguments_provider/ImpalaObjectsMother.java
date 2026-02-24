package grok_connect.providers.arguments_provider;

import grok_connect.connectors_info.FuncCall;
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
import java.time.LocalDate;
import java.time.Year;
import java.time.temporal.TemporalAdjusters;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

@SuppressWarnings("unused")
public class ImpalaObjectsMother {
    private static final Parser parser = new DateParser();

    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn("TABLE_SCHEMA", new String[] {"default"}));
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> getSchema_ok() {
        String firstColumnName = "table_schema";
        String secondColumnName = "table_name";
        String thirdColumnName = "column_name";
        String fourthColumnName = "data_type";
        String schema = "default";
        String table = "MOCK_DATA";
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn(thirdColumnName, new String[] {"id", "first_name", "last_name", "email",
                        "gender", "ip_address", "bool", "country", "dat", "some_number"}),
                new StringColumn(fourthColumnName, new String[] {"BIGINT", "STRING",
                        "STRING", "STRING", "STRING", "STRING",
                        "BOOLEAN", "STRING", "DATE", "DECIMAL(5,2)"}),
                new StringColumn(firstColumnName, new String[] {schema, schema,
                        schema, schema, schema, schema, schema,
                        schema, schema, schema}),
                new StringColumn(secondColumnName, new String[] {table, table,
                        table, table, table, table,
                        table, table, table, table}));
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> checkListParameterSupport_ok() {
        // list<string>
        List<String> values = new ArrayList<>();
        values.add("Poland");
        values.add("Brazil");
        DataFrame expected = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"2", "5", "20"}),
                new StringColumn("first_name", new String[]{"Nicholle", "Mitchell", "Lucius",}),
                new StringColumn("last_name", new String[]{"Karoly", "Haglington", "Edelmann"}),
                new StringColumn("email", new String[]{"nkaroly1@alexa.com", "mhaglington4@indiegogo.com",
                        "ledelmannj@bravesites.com"}),
                new StringColumn("gender", new String[]{"Female", "Male", "Male"}),
                new StringColumn("ip_address", new String[]{"255.233.247.118/32", "209.93.181.190/32", "66.174.30.225/32"}),
                new BoolColumn("bool", new Boolean[]{false, true, false}),
                new StringColumn("country", new String[]{"Poland", "Poland", "Brazil"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles("yyyy-MM-dd", "2014-02-27",
                                "2020-10-09","1999-06-22")),
                new FloatColumn("some_number", new Float[]{864.09f, 15.22f, 378.73f}));
        FuncCall funcCall1 = FuncCallBuilder.getBuilder()
                .addQuery("--input: list<string> values\n" +
                        "SELECT * FROM mock_data WHERE country IN (@values) ORDER BY id")
                .addFuncParam("list", "string", "values", values, "")
                .addFuncCallOptionsPattern("country", "", "",
                        null, null, "Poland", "Brazil")
                .build();
        FuncCall funcCall2 = FuncCallBuilder.getBuilder()
                .addQuery("--input: list<string> values = ['Poland','Brazil']\n" +
                        "SELECT * FROM mock_data WHERE country IN (@values) ORDER BY id")
                .addFuncParam("list", "string", "values", values, "")
                .addFuncCallOptionsPattern("country", "", "",
                        null, null, "Poland", "Brazil")
                .build();
        return Stream.of(
                Arguments.of(Named.of("type: list<string>; operator: none; pattern: none", funcCall1), expected),
                Arguments.of(Named.of("type: list<string>; operator: none; pattern: none", funcCall2), expected));
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
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "", "dat", "today", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        now.toString(), now.plusDays(1).toString())
                .build();
        // --input: string date = "this week" {pattern: datetime}
        DataFrame expected2 = DataFrame.fromColumns(
                new DateTimeColumn("dat", parser.parseDatesToDoubles(datePattern,
                                now.toString(),
                                dayOfWeek == 1 ? null : yesterday.toString(),
                                lastDayOfWeek.equals(now) ? null : lastDayOfWeek.toString())));
        FuncCall funcCall2 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"this week\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string","", "dat", "this week", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        firstDayOfWeek.toString(),
                        lastDayOfWeek.plusDays(1).toString())
                .build();
        // --input: string date = "this month" {pattern: datetime}
        DataFrame expected3 = DataFrame.fromColumns(
                new DateTimeColumn("dat", parser.parseDatesToDoubles(datePattern,
                                now.toString(),
                                dayOfMonth == 1 ? null : yesterday.toString(),
                                lastDayOfWeek.getMonthValue() >  lastDayOfMonth.getMonthValue() || lastDayOfWeek.equals(now)?
                                        null : lastDayOfWeek.toString())));
        FuncCall funcCall3 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"this month\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "","dat", "this month", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        firstDayOfMonth.toString(),
                        lastDayOfMonth.plusDays(1).toString())
                .build();
        // --input: string date = "this year" {pattern: datetime}
        DataFrame expected4 = DataFrame.fromColumns(
                new DateTimeColumn("dat", parser.parseDatesToDoubles(datePattern,
                                now.toString(),
                                dayOfMonth == 1 ? null : yesterday.toString(),
                                lastDayOfWeek.getYear() >  now.getYear() || lastDayOfWeek.equals(now)?
                                        null : lastDayOfWeek.toString(),
                                dayOfLastYear.getYear() == now.getYear() ? dayOfLastYear.toString() : null)));
        FuncCall funcCall4 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"this year\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY id\n"
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
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "","dat", "yesterday", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        yesterday.toString(),
                        now.toString())
                .build();
        // --input: string date = "last year" {pattern: datetime}
        DataFrame expected6 = DataFrame.fromColumns(
                new DateTimeColumn("dat", parser.parseDatesToDoubles(datePattern,
                                yesterday.getYear() < now.getYear() ? yesterday.toString() : null,
                                dayOfLastYear.getYear() < now.getYear() ? dayOfLastYear.toString() : null)));
        FuncCall funcCall6 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"last year\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "", "dat", "last year", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        firstDayOfYear.minusYears(1).toString(), firstDayOfYear.toString())
                .build();
        // --input: string date = "anytime" {pattern: datetime}
        DataFrame expected7 = DataFrame.fromColumns(
                new DateTimeColumn("dat", parser.parseDatesToDoubles(datePattern, now.toString(),
                                yesterday.toString(), lastDayOfWeek.getYear() >  now.getYear() || lastDayOfWeek.equals(now)?
                                        null : lastDayOfWeek.toString(), dayOfLastYear.toString(),
                                "2021-04-09")));
        FuncCall funcCall7 = FuncCallBuilder.getBuilder()
                .addQuery( "--input: string dat = \"anytime\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "", "dat", "anytime", "datetime")
                .addFuncCallOptionsPattern("dat", "", "none", true, true)
                .build();
        // --input: string date = "2021-2022" {pattern: datetime}
        DataFrame expected8 = DataFrame.fromColumns(
                new DateTimeColumn("dat", parser.parseDatesToDoubles(datePattern, "2021-04-09")));

        FuncCall funcCall8 = FuncCallBuilder.getBuilder()
                .addQuery( "--input: string dat = \"2021-2021\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "", "dat", "2021-2022", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, true,
                        Year.of(2021).atDay(1).toString(),
                        Year.of(2022).atDay(1).toString())
                .build();
        // --input: string date = "before 1/1/2022" {pattern: datetime}

        FuncCall funcCall9 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"before 1/1/2022\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "","dat", "before 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("dat", "", "before", true, true,
                        Year.of(2022).atDay(1).toString())
                .build();
        // --input: string date = "after 1/1/2022" {pattern: datetime}
        DataFrame expected9 = DataFrame.fromColumns(
                new DateTimeColumn("dat", parser.parseDatesToDoubles(datePattern, now.toString(),
                                yesterday.toString(), lastDayOfWeek.toString(), dayOfLastYear.toString())));
        FuncCall funcCall10 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"after 1/1/2022\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "","dat", "after 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("dat", "", "after", true, true,
                        LocalDate.parse("2022-01-01").toString())
                .build();
        // --input: string date = "April 2021" {pattern: datetime}
        FuncCall funcCall11 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"April 2021\" {pattern: datetime}\n"
                        + "SELECT dat FROM dates_patterns WHERE @dat(dat) ORDER BY id\n"
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

    public static Stream<Arguments> checkCharacterTypesSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn("char_type", new String[]{"Datagrok"}),
                new StringColumn("varchar_type", new String[]{"Datagrok"}),
                new StringColumn("string_type", new String[]{"Hello, world!"}));
        return Stream.of(Arguments.of(Named.of("CHARACTER TYPES SUPPORT",
                        FuncCallBuilder.fromQuery("SELECT * FROM CHARACTER_TYPES")),
                expected));
    }

    public static Stream<Arguments> checkIntegerTypesSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new BigIntColumn("bigint_type", new String[]{"9223372036854775807", "-10", "-9223372036854775808"}),
                new IntColumn("int_type", new Integer[]{2147483647, 0, -2147483648}),
                new IntColumn("smallint_type", new Integer[]{32767, -100, -32768}),
                new IntColumn("tinyint_type", new Integer[]{127, 1, -128}));
        return Stream.of(Arguments.of(Named.of("INTEGER TYPES SUPPORT",
                        FuncCallBuilder.fromQuery("SELECT * FROM INTEGER_TYPES ORDER BY tinyint_type DESC")),
                expected));
    }

    public static Stream<Arguments> checkFloatTypesSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new FloatColumn("double_type", new Float[]{Float.NEGATIVE_INFINITY, Float.NaN}),
                new FloatColumn("float_type", new Float[]{Float.POSITIVE_INFINITY, 1.401E-42f}),
                new FloatColumn("decimal_type", new Float[]{100000.0f, -1.20001E-4f}));
        return Stream.of(Arguments.of(Named.of("FLOAT TYPES SUPPORT",
                        FuncCallBuilder.fromQuery("SELECT * FROM FLOAT_TYPES ORDER BY decimal_type DESC")),
                expected));
    }

    public static Stream<Arguments> checkDateTypesSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new DateTimeColumn("date_type", parser.parseDatesToDoubles("yyyy-MM-dd", "2023-02-06",
                        "1111-11-11")),
                new DateTimeColumn("stamp_type", new Double[]{parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss",
                        "2001-01-09 01:05:01"), parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss.SSS",
                        "1985-09-25 17:45:30.005")}));
        return Stream.of(Arguments.of(Named.of("DATE TYPES SUPPORT",
                        FuncCallBuilder.fromQuery("SELECT * FROM DATE_TYPES ORDER BY date_type DESC")),
                expected));
    }
}

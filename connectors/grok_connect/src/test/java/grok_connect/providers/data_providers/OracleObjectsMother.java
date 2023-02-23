package grok_connect.providers.data_providers;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.DataFrameBuilder;
import grok_connect.providers.utils.DateParser;
import grok_connect.providers.utils.FuncCallBuilder;
import grok_connect.providers.utils.Parser;
import org.junit.jupiter.api.Named;
import org.junit.jupiter.params.provider.Arguments;
import serialization.*;

import java.time.LocalDate;
import java.time.Year;
import java.time.temporal.TemporalAdjusters;
import java.util.stream.Stream;

/**
 * Provides data for Oracle provider tests
 */
public class OracleObjectsMother {
    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[] {"DATAGROK", "GSMADMIN_INTERNAL",
                        "LBACSYS"}), "TABLE_SCHEMA")
                .build();
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> getSchema_ok() {
        String firstColumnName = "TABLE_SCHEMA";
        String secondColumnName = "TABLE_NAME";
        String thirdColumnName = "COLUMN_NAME";
        String fourthColumnName = "DATA_TYPE";
        String schema = "DATAGROK";
        String table = "MOCK_DATA";
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(9)
                .setColumn(new StringColumn(), firstColumnName, new String[] {schema, schema,
                        schema, schema, schema, schema, schema,
                        schema, schema})
                .setColumn(new StringColumn(), secondColumnName, new String[] {table, table,
                        table, table, table, table,
                        table, table, table})
                .setColumn(new StringColumn(), thirdColumnName, new String[] {"ID", "FIRST_NAME", "LAST_NAME", "EMAIL",
                        "GENDER", "IP_ADDRESS", "COUNTRY", "DAT", "SOME_NUMBER"})
                .setColumn(new StringColumn(), fourthColumnName, new String[] {"NUMBER", "VARCHAR2",
                        "VARCHAR2", "VARCHAR2", "VARCHAR2", "VARCHAR2",
                        "VARCHAR2", "DATE", "NUMBER"})
                .build();
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> checkDatesParameterSupport_ok() {
        Parser parser = new DateParser();
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
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble(datePattern, now.toString())}),
                        "dat")
                .build();
        FuncCall funcCall1 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"today\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "dat", "today", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        now.toString(), now.plusDays(1).toString())
                .build();
        // --input: string date = "this week" {pattern: datetime}
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(dayOfWeek == 1 ? 2 : 3)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern,
                                now.toString(),
                                dayOfWeek == 1 ? null : yesterday.toString(),
                                lastDayOfWeek.toString())),
                        "dat")
                .build();
        FuncCall funcCall2 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"this week\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "dat", "this week", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        firstDayOfWeek.toString(),
                        lastDayOfWeek.plusDays(1).toString())
                .build();
        // --input: string date = "this month" {pattern: datetime}
        DataFrame expected3 = DataFrameBuilder.getBuilder()
                .setRowCount(dayOfMonth > 1 && dayOfMonth < 31 - 6 ? 3 : 2)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern,
                                now.toString(),
                                dayOfMonth == 1 ? null : yesterday.toString(),
                                lastDayOfWeek.toString())),
                        "dat")
                .build();
        FuncCall funcCall3 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"this month\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "dat", "this month", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        firstDayOfMonth.toString(),
                        lastDayOfMonth.plusDays(1).toString())
                .build();
        // --input: string date = "this year" {pattern: datetime}
        DataFrame expected4 = DataFrameBuilder.getBuilder()
                .setRowCount(dayOfYear > 1 && dayOfYear < Year.of(now.getYear()).length() - 6 ? 3 : 2)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern,
                                now.toString(),
                                dayOfMonth == 1 ? null : yesterday.toString(),
                                lastDayOfWeek.toString())),
                        "dat")
                .build();
        FuncCall funcCall4 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"this year\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "dat", "this year", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        firstDayOfYear.toString(),
                        lastDayOfYear.plusDays(1).toString())
                .build();
        // --input: string date = "yesterday" {pattern: datetime}
        DataFrame expected5 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, yesterday.toString())),
                        "dat")
                .build();
        FuncCall funcCall5 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string date = \"yesterday\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "dat", "yesterday", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        yesterday.toString(),
                        now.toString())
                .build();
        // --input: string date = "last year" {pattern: datetime}
        DataFrame expected6 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern,
                                dayOfLastYear.toString())),
                        "dat")
                .build();
        FuncCall funcCall6 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"last year\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "dat", "last year", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        firstDayOfYear.minusYears(1).toString(), firstDayOfYear.toString())
                .build();
        // --input: string date = "anytime" {pattern: datetime}
        DataFrame expected7 = DataFrameBuilder.getBuilder()
                .setRowCount(5)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, now.toString(),
                                yesterday.toString(), lastDayOfWeek.toString(), dayOfLastYear.toString(),
                                "2021-04-09")),
                        "dat")
                .build();
        FuncCall funcCall7 = FuncCallBuilder.getBuilder()
                .addQuery( "--input: string dat = \"anytime\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "dat", "anytime", "datetime")
                .addFuncCallOptionsPattern("dat", "", "none", true, true)
                .build();
        // --input: string date = "2021-2022" {pattern: datetime}
        DataFrame expected8 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2021-04-09")),
                        "dat")
                .build();

        FuncCall funcCall8 = FuncCallBuilder.getBuilder()
                .addQuery( "--input: string dat = \"2021-2021\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "dat", "2021-2022", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, true,
                        Year.of(2021).atDay(1).toString(),
                        Year.of(2022).atDay(1).toString())
                .build();
        // --input: string date = "before 1/1/2022" {pattern: datetime}

        FuncCall funcCall9 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"before 1/1/2022\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "dat", "before 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("dat", "", "before", true, true,
                        Year.of(2022).atDay(1).toString())
                .build();
        // --input: string date = "after 1/1/2022" {pattern: datetime}
        DataFrame expected9 = DataFrameBuilder.getBuilder()
                .setRowCount(4)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, now.toString(),
                                yesterday.toString(), lastDayOfWeek.toString(), dayOfLastYear.toString())),
                        "dat")
                .build();
        FuncCall funcCall10 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"after 1/1/2022\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "dat", "after 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("dat", "", "after", true, true,
                        LocalDate.parse("2022-01-01").toString())
                .build();
        // --input: string date = "April 2021" {pattern: datetime}
        FuncCall funcCall11 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"April 2021\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "dat", "April 2021", "datetime")
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

    public static Stream<Arguments> checkMultipleParametersSupport() {
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
                        + "--input: string dat = \"before 1/1/2022\" {pattern: datetime}\n"
                        + "SELECT * FROM mock_data WHERE @first_name(first_name) AND @id(id) "
                        + "AND @email(email) AND @some_number(some_number) "
                        + "AND @country(country) AND @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "first_name", "starts with p", "string")
                .addFuncParam("string", "id", ">1", "string")
                .addFuncParam("string", "email", "contains com", "string")
                .addFuncParam("string", "some_number", ">20", "double")
                .addFuncParam("string", "country", "in (Indonesia)", "string")
                .addFuncParam("string", "dat", "before 1/1/2022", "datetime")
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
                .addFuncCallOptionsPattern("dat", "before 1/1/2022", "before",
                        true, true, Year.of(2022).atMonth(1).atDay(1).toString())
                .build();
        return Stream.of(Arguments.of(Named.of("type: multiple; operator: multiple; pattern: multiple", funcCall1),
                expected1));
    }
}

package grok_connect.providers.utils;

import grok_connect.connectors_info.FuncCall;
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
import java.time.LocalDate;
import java.time.Year;
import java.time.temporal.TemporalAdjusters;
import java.util.Base64;
import java.util.stream.Stream;

/**
 * Provides data for tests
 */
public class ObjectsMother {
    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[] {"public", "pg_catalog",
                        "information_schema"}), "table_schema")
                .build();
        return Stream.of(Arguments.of(expected));
    }
    public static Stream<Arguments> checkParameterSupport_ok() {
        Parser parser = new DateParser();
        String datePattern = "yyyy-MM-dd";
        // --input: int id = 20
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new BigIntColumn(new String[]{"20"}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Lucius"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Edelmann"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"ledelmannj@bravesites.com"}), "email")
                .setColumn(new StringColumn(new String[]{"Male"}), "gender")
                .setColumn(new StringColumn(new String[]{"66.174.30.225/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{false}), "bool")
                .setColumn(new StringColumn(new String[]{"Brazil"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "1999-06-22")), "date")
                .setColumn(new FloatColumn(new Float[]{378.73f}), "some_number")
                .build();
        FuncCall funcCall1 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByInt\n"
                        + "--input: int id = 20\n"
                        + "SELECT * FROM mock_data WHERE id = @id;\n"
                        + "--end")
                .addFuncParam("int", "id", 20, "")
                .build();
        // --input: string id = ">28" {pattern: int}
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new BigIntColumn(new String[]{"29", "30"}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Grantham", "Bran"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Fayter", "Longlands"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"gfayters@desdev.cn", "blonglandst@tripod.com"}), "email")
                .setColumn(new StringColumn(new String[]{"Male", "Genderqueer"}), "gender")
                .setColumn(new StringColumn(new String[]{"26.120.76.78/32", "14.92.3.30/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{false, false}), "bool")
                .setColumn(new StringColumn(new String[]{"Sweden", "France"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2009-10-02",
                        "2016-07-10")), "date")
                .setColumn(new FloatColumn(new Float[]{595.22f, 879.94f}), "some_number")
                .build();

        FuncCall funcCall2 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternInt\n"
                        + "--input: string id = \">28\" {pattern: int}\n"
                        + "SELECT * FROM mock_data WHERE @id(id)\n"
                        + "--end")
                .addFuncParam("string", "id", ">28", "int")
                .addFuncCallOptionsPattern("id", ">28", ">",
                        null, null, 28)
                .build();
        // input: string id = ">=29" {pattern: int}
        FuncCall funcCall3 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternInt\n"
                        + "--input: string id = \">=29\" {pattern: int}\n"
                        + "SELECT * FROM mock_data WHERE @id(id)\n"
                        + "--end")
                .addFuncParam("string", "id", ">=29", "int")
                .addFuncCallOptionsPattern("id", ">=29", ">=",
                        null, null, 29)
                .build();
        // --input: string id = "<=1" {pattern: int}
        DataFrame expected3 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new BigIntColumn(new String[]{"1"}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Burk"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Kemery"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"bkemery0@businesswire.com"}), "email")
                .setColumn(new StringColumn(new String[]{"Male"}), "gender")
                .setColumn(new StringColumn(new String[]{"249.64.22.121/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{true}), "bool")
                .setColumn(new StringColumn(new String[]{"China"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2017-09-20")), "date")
                .setColumn(new FloatColumn(new Float[]{510.32f}), "some_number")
                .build();
        FuncCall funcCall4 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternInt\n"
                        + "--input: string id = \"<=1\" {pattern: int}\n"
                        + "SELECT * FROM mock_data WHERE @id(id)\n"
                        + "--end")
                .addFuncParam("string", "id", "<=1", "int")
                .addFuncCallOptionsPattern("id", "<=1", "<=",
                        null, null, 1)
                .build();
        // --input: string id = "<2" {pattern: int}
        FuncCall funcCall5 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternInt\n"
                        + "--input: string id = \"<2\" {pattern: int}\n"
                        + "SELECT * FROM mock_data WHERE @id(id)\n"
                        + "--end")
                .addFuncParam("string", "id", "<2", "int")
                .addFuncCallOptionsPattern("id", "<2", "<",
                        null, null, 2)
                .build();
        // --input: string id = "in(29, 30)" {pattern: int}
        FuncCall funcCall6 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternInt\n"
                        + "--input: string id = \"in(29, 30)\" {pattern: int}\n"
                        + "SELECT * FROM mock_data WHERE @id(id)\n"
                        + "--end")
                .addFuncParam("string", "id", "in(29, 30)", "int")
                .addFuncCallOptionsPattern("id", "in(29, 30)", "in",
                        null, null, 29, 30)
                .build();
        // --input: string id = "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
        DataFrame expected4 = DataFrameBuilder.getBuilder()
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
        FuncCall funcCall7 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternInt\n"
                        + "--input: string id = \"not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)\" {pattern: int}\n"
                        + "SELECT * FROM mock_data WHERE @id(id)\n"
                        + "--end")
                .addFuncParam("string", "id", "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)",
                        "int")
                .addFuncCallOptionsPattern("id", "not in(21, 22, 23, 24, 25, 26, 27, 28, 29, 30)",
                        "not in", null, null, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)
                .build();
        // --input: string id = "min-max 29-30" {pattern: int}
        FuncCall funcCall8 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternInt\n"
                        + "--input: string id = \"min-max 29-30\" {pattern: int}\n"
                        + "SELECT * FROM mock_data WHERE @id(id)\n"
                        + "--end")
                .addFuncParam("string", "id", "min-max 29-30",
                        "int")
                .addFuncCallOptionsPattern("id", "29-30",
                        "-", null, null, 29, 30)
                .build();
        //--input: double some_number = 510.32
        FuncCall funcCall9 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByDouble\n"
                        + "--input: double some_number = 510.32\n"
                        + "SELECT * FROM mock_data WHERE some_number = @some_number;\n"
                        + "--end")
                .addFuncParam("double", "some_number", 510.32, "double")
                .build();
        // --input: string some_number = ">975" {pattern: double}
        DataFrame expected5 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new BigIntColumn(new String[]{"10", "26"}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Scottie", "Daryle"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Formilli", "O'Shaughnessy"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"sformilli9@aol.com", "doshaughnessyp@com.com"}),
                        "email")
                .setColumn(new StringColumn(new String[]{"Male", "Male"}), "gender")
                .setColumn(new StringColumn(new String[]{"101.241.191.228/32", "204.107.16.207/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{false, false}), "bool")
                .setColumn(new StringColumn(new String[]{"Vietnam", "Honduras"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2003-01-04",
                        "2010-05-04")), "date")
                .setColumn(new FloatColumn(new Float[]{978.01f, 983.03f}), "some_number")
                .build();
        FuncCall funcCall10 =  FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByDouble\n" +
                        "--input: string some_number = \">975\" {pattern: double}\n"
                        + "SELECT * FROM mock_data WHERE @some_number(some_number);\n"
                        + "--end")
                .addFuncParam("string", "some_number", ">975", "double")
                .addFuncCallOptionsPattern("some_number", ">975", ">",
                        null, null, 975)
                .build();
        // --input: string some_number = ">=975" {pattern: double}
        FuncCall funcCall11 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByDouble\n"
                        + "--input: string some_number = \">=975\" {pattern: double}\n"
                        + "SELECT * FROM mock_data WHERE @some_number(some_number);\n"
                        + "--end")
                .addFuncParam("string", "some_number", ">=975", "double")
                .addFuncCallOptionsPattern("some_number", ">=975", ">=",
                        null, null, 975)
                .build();
        //--input: string some_number = "<20" {pattern: double}
        DataFrame expected6 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new BigIntColumn(new String[]{"5"}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Mitchell"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Haglington"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"mhaglington4@indiegogo.com"}), "email")
                .setColumn(new StringColumn(new String[]{"Male"}), "gender")
                .setColumn(new StringColumn(new String[]{"209.93.181.190/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{true}), "bool")
                .setColumn(new StringColumn(new String[]{"Poland"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2020-10-09")), "date")
                .setColumn(new FloatColumn(new Float[]{15.22f}), "some_number")
                .build();
        FuncCall funcCall12 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByDouble\n"
                        + "--input: string some_number = \"<20\" {pattern: double}\n"
                        + "SELECT * FROM mock_data WHERE @some_number(some_number);\n"
                        + "--end")
                .addFuncParam("string", "some_number", "<20", "double")
                .addFuncCallOptionsPattern("some_number", "<20", "<",
                        null, null, 20)
                .build();
        // --input: string some_number = "<=20" {pattern: double}
        FuncCall funcCall13 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByDouble\n" +
                        "--input: string some_number = \"<=20\" {pattern: double}\n"
                        + "SELECT * FROM mock_data WHERE @some_number(some_number);\n"
                        + "--end")
                .addFuncParam("string", "some_number", "<=20", "double")
                .addFuncCallOptionsPattern("some_number", "<=20", "<=",
                        null, null, 20)
                .build();
        // --input: string first_name = 'contains Z' {pattern: string}
        DataFrame expected7 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new BigIntColumn(new String[]{"25"}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Zolly"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Wimmers"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"zwimmerso@hatena.ne.jp"}), "email")
                .setColumn(new StringColumn(new String[]{"Male"}), "gender")
                .setColumn(new StringColumn(new String[]{"123.12.225.114/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{false}), "bool")
                .setColumn(new StringColumn(new String[]{"Bosnia and Herzegovina"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2003-02-12")), "date")
                .setColumn(new FloatColumn(new Float[]{217.18f}), "some_number")
                .build();
        FuncCall funcCall14 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternString\n"
                        + "--input: string first_name = 'contains Z' {pattern: string}\n"
                        + "SELECT * FROM mock_data WHERE @first_name(first_name);\n"
                        + "--end")
                .addFuncParam("string", "first_name", "contains Z", "string")
                .addFuncCallOptionsPattern("first_name", "contains Z", "contains",
                        null, null, "Z")
                .build();
        // --input: string first_name = 'starts with W' {pattern: string}
        DataFrame expected8 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new BigIntColumn(new String[]{"23"}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Waly"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Rogliero"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"wroglierom@berkeley.edu"}), "email")
                .setColumn(new StringColumn(new String[]{"Female"}), "gender")
                .setColumn(new StringColumn(new String[]{"122.90.196.231/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{true}), "bool")
                .setColumn(new StringColumn(new String[]{"Sweden"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2011-12-18")),
                        "date")
                .setColumn(new FloatColumn(new Float[]{147.69f}), "some_number")
                .build();
        FuncCall funcCall15 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternString\n"
                        + "--input: string first_name = 'starts with W' {pattern: string}\n"
                        + "SELECT * FROM mock_data WHERE @first_name(first_name);\n"
                        + "--end")
                .addFuncParam("string", "first_name", "starts with W", "string")
                .addFuncCallOptionsPattern("first_name", "starts with W", "starts with",
                        null, null, "W")
                .build();
        // --input: string first_name = 'ends with y' {pattern: string}
        DataFrame expected9 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new BigIntColumn(new String[]{"20"}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Lucius"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Edelmann"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"ledelmannj@bravesites.com"}), "email")
                .setColumn(new StringColumn(new String[]{"Male"}), "gender")
                .setColumn(new StringColumn(new String[]{"66.174.30.225/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{false}), "bool")
                .setColumn(new StringColumn(new String[]{"Brazil"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "1999-06-22")),
                        "date")
                .setColumn(new FloatColumn(new Float[]{378.73f}), "some_number")
                .build();
        FuncCall funcCall16 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternString\n"
                        + "--input: string first_name = 'ends with s' {pattern: string}\n"
                        + "SELECT * FROM mock_data WHERE @first_name(first_name);\n"
                        + "--end")
                .addFuncParam("string", "first_name", "ends with s", "string")
                .addFuncCallOptionsPattern("first_name", "ends with s", "ends with",
                        null, null, "s")
                .build();
        // --input: string country = 'in (Poland, Brazil)' {pattern: string}
        DataFrame expected10 = DataFrameBuilder.getBuilder()
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
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2014-02-27",
                                "2020-10-09","1999-06-22")),
                        "date")
                .setColumn(new FloatColumn(new Float[]{864.09f, 15.22f, 378.73f}), "some_number")
                .build();
        FuncCall funcCall17 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternString\n" +
                        "--input: string country = 'in (Poland, Brazil)' {pattern: string}\n" +
                        "SELECT * FROM mock_data WHERE @country(country);\n" +
                        "--end")
                .addFuncParam("string", "country", "in (Poland, Brazil)", "string")
                .addFuncCallOptionsPattern("country", "in (Poland, Brazil)", "in",
                        null, null, "Poland", "Brazil")
                .build();
        // --input: string email = 'regex ^([A-Za-z0-9_]+@google.com.au)$' {pattern: string}
        DataFrame expected11 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new BigIntColumn(new String[]{"9"}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Marlie"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Mayze"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"mmayze8@google.com.au"}), "email")
                .setColumn(new StringColumn(new String[]{"Female"}), "gender")
                .setColumn(new StringColumn(new String[]{"68.41.25.65/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{false}), "bool")
                .setColumn(new StringColumn(new String[]{"France"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2011-11-10")),
                        "date")
                .setColumn(new FloatColumn(new Float[]{561.72f}), "some_number")
                .build();
        FuncCall funcCall18 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternString\n"
                        + "--input: string email = 'regex ^([A-Za-z0-9_]+@google.com.au)$' {pattern: string}\n"
                        + "SELECT * FROM mock_data WHERE @email(email);\n"
                        + "--end")
                .addFuncParam("string", "email", "regex ^([A-Za-z0-9_]+@google.com.au)$", "string")
                .addFuncCallOptionsPattern("email", "regex ^([A-Za-z0-9_]+@google.com.au)$",
                        "regex", null, null, "^([A-Za-z0-9_]+@google.com.au)$")
                .build();
        // --input: string first_name = "starts with p" {pattern: string}
        //--input: string id = ">1" {pattern :int}
        //--input: bool bool = false
        //--input: string email = "contains com" {pattern: string}
        //--input: string some_number = ">20" {pattern: double}
        //--input: string country = "in (Indonesia)" {pattern: string}
        //--input: string date = "before 1/1/2022" {pattern: datetime}
        DataFrame expected12 = DataFrameBuilder.getBuilder()
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
        FuncCall funcCall19 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlAll\n" +
                        "--input: string first_name = \"starts with p\" {pattern: string}\n" +
                        "--input: string id = \">1\" {pattern :int}\n" +
                        "--input: bool bool = false\n" +
                        "--input: string email = \"contains com\" {pattern: string}\n" +
                        "--input: string some_number = \">20\" {pattern: double}\n" +
                        "--input: string country = \"in (Indonesia)\" {pattern: string}\n" +
                        "--input: string date = \"before 1/1/2022\" {pattern: datetime}\n" +
                        "SELECT * FROM mock_data\n" +
                        "WHERE @first_name(first_name)\n" +
                        "  AND @id(id)\n" +
                        "           AND bool = @bool\n" +
                        "           AND @email(email)\n" +
                        "           AND @some_number(some_number)\n" +
                        "           AND @country(country)\n" +
                        "           AND @date(date);\n" +
                        "\n" +
                        "--end")
                .addFuncParam("string", "first_name", "starts with p", "string")
                .addFuncParam("string", "id", ">1", "string")
                .addFuncParam("bool", "bool", false, "")
                .addFuncParam("string", "email", "contains com", "string")
                .addFuncParam("string", "some_number", ">20", "double")
                .addFuncParam("string", "country", "in (Indonesia)", "string")
                .addFuncParam("string", "date", "before 1/1/2022", "datetime")
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
        return Stream.of(
                Arguments.of(Named.of("type: int; operator: =; pattern: none", funcCall1), expected1),
                Arguments.of(Named.of("type: string; operator: >; pattern: int", funcCall2), expected2),
                Arguments.of(Named.of("type: string; operator: >=; pattern: int", funcCall3), expected2),
                Arguments.of(Named.of("type: string; operator: <=; pattern: int", funcCall4), expected3),
                Arguments.of(Named.of("type: string; operator: <; pattern: int", funcCall5), expected3),
                Arguments.of(Named.of("type: string; operator: in; pattern: int", funcCall6), expected2),
                Arguments.of(Named.of("type: string; operator: not in; pattern: int", funcCall7), expected4),
                Arguments.of(Named.of("type: string; operator: min-max; pattern: int", funcCall8), expected2),
                Arguments.of(Named.of("type: double; operator: =; pattern: none", funcCall9), expected3),
                Arguments.of(Named.of("type: string; operator: >; pattern: double", funcCall10), expected5),
                Arguments.of(Named.of("type: string; operator: >=; pattern: double", funcCall11), expected5),
                Arguments.of(Named.of("type: string; operator: <; pattern: double", funcCall12), expected6),
                Arguments.of(Named.of("type: string; operator: <=; pattern: double", funcCall13), expected6),
                Arguments.of(Named.of("type: string; operator: contains; pattern: string", funcCall14), expected7),
                Arguments.of(Named.of("type: string; operator: starts with; pattern: string", funcCall15), expected8),
                Arguments.of(Named.of("type: string; operator: ends with; pattern: string", funcCall16), expected9),
                Arguments.of(Named.of("type: string; operator: in; pattern: string", funcCall17), expected10),
                Arguments.of(Named.of("type: string; operator: regex; pattern: string", funcCall18), expected11),
                Arguments.of(Named.of("type: multiple; operator: multiple; pattern: multiple", funcCall19),
                        expected12)
        );
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
                        "date")
                .build();
        FuncCall funcCall1 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternDatetime\n"
                        + "--input: string date = \"today\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date);\n"
                        + "--end")
                .addFuncParam("string", "date", "today", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        now.toString(), now.plusDays(1).toString())
                .build();
        // --input: string date = "this week" {pattern: datetime}
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(dayOfWeek == 1 ? 2 : 3)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern,
                                now.toString(),
                                dayOfWeek == 1 ? null : yesterday.toString(),
                                lastDayOfWeek.toString())),
                        "date")
                .build();
        FuncCall funcCall2 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternDatetime\n"
                        + "--input: string date = \"this week\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date);\n"
                        + "--end")
                .addFuncParam("string", "date", "this week", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
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
                        "date")
                .build();
        FuncCall funcCall3 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternDatetime\n"
                        + "--input: string date = \"this month\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date);\n"
                        + "--end")
                .addFuncParam("string", "date", "this month", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
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
                        "date")
                .build();
        FuncCall funcCall4 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternDatetime\n"
                        + "--input: string date = \"this year\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date);\n"
                        + "--end")
                .addFuncParam("string", "date", "this year", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        firstDayOfYear.toString(),
                        lastDayOfYear.plusDays(1).toString())
                .build();
        // --input: string date = "yesterday" {pattern: datetime}
        DataFrame expected5 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, yesterday.toString())),
                        "date")
                .build();
        FuncCall funcCall5 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternDatetime\n"
                        + "--input: string date = \"yesterday\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date);\n"
                        + "--end")
                .addFuncParam("string", "date", "yesterday", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        yesterday.toString(),
                        now.toString())
                .build();
        // --input: string date = "last year" {pattern: datetime}
        DataFrame expected6 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern,
                               dayOfLastYear.toString())),
                        "date")
                .build();
        FuncCall funcCall6 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternDatetime\n"
                        + "--input: string date = \"last year\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date);\n"
                        + "--end")
                .addFuncParam("string", "date", "last year", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        firstDayOfYear.minusYears(1).toString(), firstDayOfYear.toString())
                .build();
        // --input: string date = "anytime" {pattern: datetime}
        DataFrame expected7 = DataFrameBuilder.getBuilder()
                .setRowCount(5)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, now.toString(),
                                yesterday.toString(), lastDayOfWeek.toString(), dayOfLastYear.toString(),
                                "2021-04-09")),
                        "date")
                .build();
        FuncCall funcCall7 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternDatetime\n"
                        + "--input: string date = \"anytime\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date);\n"
                        + "--end")
                .addFuncParam("string", "date", "anytime", "datetime")
                .addFuncCallOptionsPattern("date", "", "none", true, true)
                .build();
        // --input: string date = "2021-2022" {pattern: datetime}
        DataFrame expected8 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2021-04-09")),
                        "date")
                .build();

        FuncCall funcCall8 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternDatetime\n"
                        + "--input: string date = \"2021-2021\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date);\n"
                        + "--end")
                .addFuncParam("string", "date", "2021-2022", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, true,
                        Year.of(2021).atDay(1).toString(),
                        Year.of(2022).atDay(1).toString())
                .build();
        // --input: string date = "before 1/1/2022" {pattern: datetime}

        FuncCall funcCall9 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternDatetime\n"
                        + "--input: string date = \"before 1/1/2022\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date);\n"
                        + "--end")
                .addFuncParam("string", "date", "before 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("date", "", "before", true, true,
                        Year.of(2022).atDay(1).toString())
                .build();
        // --input: string date = "after 1/1/2022" {pattern: datetime}
        DataFrame expected9 = DataFrameBuilder.getBuilder()
                .setRowCount(4)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, now.toString(),
                                yesterday.toString(), lastDayOfWeek.toString(), dayOfLastYear.toString())),
                        "date")
                .build();
        FuncCall funcCall10 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternDatetime\n"
                        + "--input: string date = \"after 1/1/2022\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date);\n"
                        + "--end")
                .addFuncParam("string", "date", "after 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("date", "", "after", true, true,
                        LocalDate.parse("2022-01-01").toString())
                .build();
        // --input: string date = "April 2021" {pattern: datetime}
        FuncCall funcCall11 = FuncCallBuilder.getBuilder()
                .addQuery("--name: PostgresqlByStringPatternDatetime\n"
                        + "--input: string date = \"April 2021\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date);\n"
                        + "--end")
                .addFuncParam("string", "date", "April 2021", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
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

    public static Stream<Arguments> getSchema_ok() {
        // based on test script postgres_basic_types.sql
        String firstColumnName = "table_schema";
        String secondColumnName = "table_name";
        String thirdColumnName = "column_name";
        String fourthColumnName = "data_type";
        String schema = "public";
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
                .setColumn(new StringColumn(), fourthColumnName, new String[] {"bigint", "character varying",
                        "character varying", "character varying", "character varying", "cidr",
                        "boolean", "character varying", "date", "numeric"})
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
                        getShortFuncCall("SELECT * FROM sal_emp;")), expected1),
                Arguments.of(Named.of("ARRAY'S INTERNAL TYPE SUPPORT#1",
                        getShortFuncCall("SELECT pay_by_quarter[3] FROM sal_emp;")), expected2),
                Arguments.of(Named.of("ARRAY'S INTERNAL TYPE SUPPORT#2",
                        getShortFuncCall("SELECT schedule[:2][2:] FROM sal_emp WHERE "
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
                getShortFuncCall("SELECT * FROM mock_data LIMIT 10;")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_bitStringType_ok() {
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new StringColumn(new String[] {"101", "001"}), "a")
                .setColumn(new StringColumn(new String[] {"0011", "101"}), "b")
                .build();
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new BoolColumn(new Boolean[] {true, false}), "c")
                .setColumn(new StringColumn(new String[] {"101", "0"}), "d")
                .build();
        return Stream.of(Arguments.of(Named.of("BIT STRING WITH FIXED LENGTH SUPPORT",
                        getShortFuncCall("SELECT * FROM test1;")), expected1),
                Arguments.of(Named.of("BIT STRING WITH VARYING LENGTH SUPPORT",
                        getShortFuncCall("SELECT * FROM test2")), expected2));
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
                getShortFuncCall("SELECT * FROM bytea_data;")), expected));
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
                        getShortFuncCall("SELECT * FROM on_hand;")), expected1),
                Arguments.of(Named.of("COMPOSITE TYPE FIELDS TYPES SUPPORT",
                        getShortFuncCall("SELECT (item).name, (item).supplier_id, "
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
                getShortFuncCall("SELECT * FROM dates")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_jsonbType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(
                        new String[]{"{\"phones\": [{\"type\": \"mobile\", \"phone\": \"001001\"}, {\"type\": \"fix\", \"phone\": \"002002\"}]}",
                                "{\"bar\": \"baz\", \"active\": false, \"balance\": 7.77}", "{\"reading\": 0.00001230}"}), "data")
                .build();
        return Stream.of(Arguments.of(Named.of("JSONB TYPE SUPPORT",
                getShortFuncCall("SELECT * FROM jsonb_data")), expected));
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
                .setColumn(new FloatColumn(new Float[] {2.2221f, 2222.2f}), "data")
                .build();
        DataFrame expected6 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new FloatColumn(new Float[] {22.22f, 55.55f}), "data")
                .build();
        return Stream.of(
                Arguments.of(Named.of("NUMERIC TYPE SUPPORT WITHOUT PRECISION AND SCALE",
                        getShortFuncCall("SELECT * FROM numeric_data")), expected1),
                Arguments.of(Named.of("DOUBLE TYPE SUPPORT",
                        getShortFuncCall("SELECT * FROM doubles")), expected2),
                Arguments.of(Named.of("REAL TYPE SUPPORT",
                        getShortFuncCall("SELECT * FROM reals")), expected3),
                Arguments.of(Named.of("BIGINT TYPE SUPPORT",
                        getShortFuncCall("SELECT * FROM bigint_data")), expected4),
                Arguments.of(Named.of("NUMERIC TYPE SUPPORT WITH PRECISION",
                        getShortFuncCall("SELECT * FROM numeric_data_precision")), expected5),
                Arguments.of(Named.of("NUMERIC TYPE SUPPORT WITH PRECISION AND SCALE",
                        getShortFuncCall("SELECT * FROM numeric_data_precision_scale")), expected6)
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
                Arguments.of(Named.of("SMALL SERIAL TYPE SUPPORT", getShortFuncCall("SELECT * FROM cars_small")), expected1),
                Arguments.of(Named.of("SERIAL TYPE SUPPORT", getShortFuncCall("SELECT * FROM cars")), expected1),
                Arguments.of(Named.of("BIG SERIAL TYPE SUPPORT", getShortFuncCall("SELECT * FROM cars_big")), expected2)
        );
    }

    public static Stream<Arguments> checkOutputDataFrame_uuidType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[]{"6ecd8c99-4036-403d-bf84-cf8400f67836",
                        "c81d4e2e-bcf2-11e6-869b-7df92533d2db", "237e9877-e79b-12d4-a765-321741963000"}), "id")
                .build();
        return Stream.of(Arguments.of(Named.of("UUID TYPE SUPPORT", getShortFuncCall("SELECT * FROM uuid_data")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_xmlType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[]{"<foo>Hello World!</foo>",
                        "abc<foo>bar</foo><bar>foo</bar",
                        "<?xml version=\"1.0\"?><book><title>Manual</title><chapter>...</chapter></book>"}), "data")
                .build();
        return Stream.of(Arguments.of(
                Named.of("XML TYPE SUPPORT", getShortFuncCall("SELECT * FROM xml_data")), expected
            )
        );
    }

    private static FuncCall getShortFuncCall(String sqlQuery) {
        return FuncCallBuilder.getBuilder()
                .addQuery(sqlQuery)
                .build();
    }
}

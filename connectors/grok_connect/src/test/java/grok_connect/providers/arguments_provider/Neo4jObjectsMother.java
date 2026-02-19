package grok_connect.providers.arguments_provider;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.DateParser;
import grok_connect.providers.utils.FuncCallBuilder;
import grok_connect.providers.utils.Parser;
import org.junit.jupiter.api.Named;
import org.junit.jupiter.params.provider.Arguments;
import serialization.*;

import java.time.LocalDate;
import java.time.Year;
import java.time.temporal.TemporalAdjusters;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

@SuppressWarnings("unused")
public class Neo4jObjectsMother {
    public static Parser parser = new DateParser();

    public static Stream<Arguments> checkParameterSupport_ok() {
        String datePattern = "yyyy-MM-dd";
        // --input: int id = 20
        DataFrame expected1 = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"20"}),
                new StringColumn("first_name", new String[]{"Lucius"}),
                new StringColumn("last_name", new String[]{"Edelmann"}),
                new StringColumn("email", new String[]{"ledelmannj@bravesites.com"}),
                new StringColumn("gender", new String[]{"Male"}),
                new StringColumn("ip_address", new String[]{"66.174.30.225/32"}),
                new BoolColumn("bool", new Boolean[]{false}),
                new StringColumn("country", new String[]{"Brazil"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "1999-06-22")),
                new FloatColumn("some_number", new Float[]{378.73f}));
        FuncCall funcCall1 = FuncCallBuilder.getBuilder()
                .addQuery("// input: int id = 20\n"
                        + "MATCH(p:Person) WHERE p.id = @id RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("int", "", "id", 20, "")
                .build();
        // --input: string id = ">28" {pattern: int}
        DataFrame expected2 = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"29", "30"}),
                new StringColumn("first_name", new String[]{"Grantham", "Bran"}),
                new StringColumn("last_name", new String[]{"Fayter", "Longlands"}),
                new StringColumn("email", new String[]{"gfayters@desdev.cn", "blonglandst@tripod.com"}),
                new StringColumn("gender", new String[]{"Male", "Genderqueer"}),
                new StringColumn("ip_address", new String[]{"26.120.76.78/32", "14.92.3.30/32"}),
                new BoolColumn("bool", new Boolean[]{false, false}),
                new StringColumn("country", new String[]{"Sweden", "France"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "2009-10-02",
                        "2016-07-10")),
                new FloatColumn("some_number", new Float[]{595.22f, 879.94f}));

        FuncCall funcCall2 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string id = \">28\" {pattern: int}\n"
                        + "MATCH(p:Person) WHERE @id(p.id) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string", "", "id", ">28", "int")
                .addFuncCallOptionsPattern("id", ">28", ">",
                        null, null, 28)
                .build();
        // input: string id = ">=29" {pattern: int}
        FuncCall funcCall3 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string id = \">=29\" {pattern: int}\n"
                        + "MATCH(p:Person) WHERE @id(p.id) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string", "", "id", ">=29", "int")
                .addFuncCallOptionsPattern("id", ">=29", ">=",
                        null, null, 29)
                .build();
        // --input: string id = "<=1" {pattern: int}
        DataFrame expected3 = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"1"}),
                new StringColumn("first_name", new String[]{"Burk"}),
                new StringColumn("last_name", new String[]{"Kemery"}),
                new StringColumn("email", new String[]{"bkemery0@businesswire.com"}),
                new StringColumn("gender", new String[]{"Male"}),
                new StringColumn("ip_address", new String[]{"249.64.22.121/32"}),
                new BoolColumn("bool", new Boolean[]{true}),
                new StringColumn("country", new String[]{"China"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "2017-09-20")),
                new FloatColumn("some_number", new Float[]{510.32f}));
        FuncCall funcCall4 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string id = \"<=1\" {pattern: int}\n"
                        + "MATCH(p:Person) WHERE @id(p.id) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string", "", "id", "<=1", "int")
                .addFuncCallOptionsPattern("id", "<=1", "<=",
                        null, null, 1)
                .build();
        // --input: string id = "<2" {pattern: int}
        FuncCall funcCall5 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string id = \"<2\" {pattern: int}\n"
                        + "MATCH(p:Person) WHERE @id(p.id) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string", "", "id", "<2", "int")
                .addFuncCallOptionsPattern("id", "<2", "<",
                        null, null, 2)
                .build();
        // --input: string id = "in(29, 30)" {pattern: int}
        FuncCall funcCall6 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string id = \"in(29, 30)\" {pattern: int}\n"
                        + "MATCH(p:Person) WHERE @id(p.id) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string", "","id", "in(29, 30)", "int")
                .addFuncCallOptionsPattern("id", "in(29, 30)", "in",
                        null, null, 29, 30)
                .build();
        // --input: string id = "not in(11, 12, 13, 14, 15, "
        //                                + "16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
        DataFrame expected4 = DataFrame.fromColumns(
                new BigIntColumn("id",
                                Stream.iterate(1, x -> x + 1).limit(10)
                                        .map(String::valueOf).toArray(String[]::new)),
                new StringColumn("first_name", new String[]{"Burk", "Nicholle", "Orlando", "Gothart",
                        "Mitchell", "Jeromy", "Joela", "Darren", "Marlie", "Scottie"}),
                new StringColumn("last_name", new String[]{"Kemery", "Karoly", "Westgate",
                                "Cokayne", "Haglington", "Twinn", "Cornau", "Juares", "Mayze", "Formilli"}),
                new StringColumn("email", new String[]{"bkemery0@businesswire.com", "nkaroly1@alexa.com",
                        "owestgate2@dedecms.com", "gcokayne3@plala.or.jp", "mhaglington4@indiegogo.com",
                        "jtwinn5@globo.com", "jcornau6@imgur.com", "djuares7@hexun.com",
                        "mmayze8@google.com.au", "sformilli9@aol.com"}),
                new StringColumn("gender", new String[]{"Male", "Female", "Polygender", "Male", "Male", "Male",
                        "Female", "Male", "Female", "Male"}),
                new StringColumn("ip_address", new String[]{"249.64.22.121/32", "255.233.247.118/32",
                                "75.0.252.254/32", "196.83.12.163/32", "209.93.181.190/32", "25.13.2.132/32",
                                "195.47.88.236/32", "94.170.16.96/32", "68.41.25.65/32", "101.241.191.228/32"}),
                new BoolColumn("bool", new Boolean[]{true, false, false, true, true, true,
                        false, false, false, false}),
                new StringColumn("country", new String[]{"China", "Poland", "Netherlands",
                        "Philippines", "Poland", "Serbia", "Indonesia", "China",
                        "France", "Vietnam"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern,"2017-09-20",
                        "2014-02-27", "2020-09-03", "2001-01-31", "2020-10-09",
                        "2014-10-04", "2020-03-19", "2011-04-09", "2011-11-10", "2003-01-04")),
                new FloatColumn("some_number", new Float[]{510.32f, 864.09f, 822.7f, 251.05f, 15.22f,
                        378.4f, 349.11f, 631.89f, 561.72f, 978.01f}));
        FuncCall funcCall7 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string id = \"not in(11, 12, 13, 14, 15, 16, 17, 18, 19, 20, "
                        + "21, 22, 23, 24, 25, 26, 27, 28, 29, 30)\" {pattern: int}\n"
                        + "MATCH(p:Person) WHERE @id(p.id) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string", "", "id", "not in(11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, " +
                                "22, 23, 24, 25, 26, 27, 28, 29, 30)",
                        "int")
                .addFuncCallOptionsPattern("id", "not in(11, 12, 13, 14, 15, "
                                + "16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)",
                        "not in", null, null, 11, 12, 13, 14, 15, 16, 17,
                        18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)
                .build();
        // --input: string id = "min-max 29-30" {pattern: int}
        FuncCall funcCall8 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string id = \"min-max 29-30\" {pattern: int}\n"
                        + "MATCH(p:Person) WHERE @id(p.id) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string", "", "id", "min-max 29-30",
                        "int")
                .addFuncCallOptionsPattern("id", "29-30",
                        "-", null, null, 29, 30)
                .build();
        //--input: double some_number = 510.32
        FuncCall funcCall9 = FuncCallBuilder.getBuilder()
                .addQuery("// input: double some_number = 510.32\n"
                        + "MATCH(p:Person) WHERE p.some_number = @some_number RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("double", "", "some_number", 510.32, "double")
                .build();
        // --input: string some_number = ">975" {pattern: double}
        DataFrame expected5 = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"10", "26"}),
                new StringColumn("first_name", new String[]{"Scottie", "Daryle"}),
                new StringColumn("last_name", new String[]{"Formilli", "O'Shaughnessy"}),
                new StringColumn("email", new String[]{"sformilli9@aol.com", "doshaughnessyp@com.com"}),
                new StringColumn("gender", new String[]{"Male", "Male"}),
                new StringColumn("ip_address", new String[]{"101.241.191.228/32", "204.107.16.207/32"}),
                new BoolColumn("bool", new Boolean[]{false, false}),
                new StringColumn("country", new String[]{"Vietnam", "Honduras"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "2003-01-04",
                        "2010-05-04")),
                new FloatColumn("some_number", new Float[]{978.01f, 983.03f}));
        FuncCall funcCall10 =  FuncCallBuilder.getBuilder()
                .addQuery("// input: string some_number = \">975\" {pattern: double}\n"
                        + "MATCH(p:Person) WHERE @some_number(p.some_number) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string",  "", "some_number", ">975", "double")
                .addFuncCallOptionsPattern("some_number", ">975", ">",
                        null, null, 975)
                .build();
        // --input: string some_number = ">=975" {pattern: double}
        FuncCall funcCall11 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string some_number = \">=975\" {pattern: double}\n"
                        + "MATCH(p:Person) WHERE @some_number(p.some_number) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string", "", "some_number", ">=975", "double")
                .addFuncCallOptionsPattern("some_number", ">=975", ">=",
                        null, null, 975)
                .build();
        //--input: string some_number = "<20" {pattern: double}
        DataFrame expected6 = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"5"}),
                new StringColumn("first_name", new String[]{"Mitchell"}),
                new StringColumn("last_name", new String[]{"Haglington"}),
                new StringColumn("email", new String[]{"mhaglington4@indiegogo.com"}),
                new StringColumn("gender", new String[]{"Male"}),
                new StringColumn("ip_address", new String[]{"209.93.181.190/32"}),
                new BoolColumn("bool", new Boolean[]{true}),
                new StringColumn("country", new String[]{"Poland"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "2020-10-09")),
                new FloatColumn("some_number", new Float[]{15.22f}));
        FuncCall funcCall12 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string some_number = \"<20\" {pattern: double}\n"
                        + "MATCH(p:Person) WHERE @some_number(p.some_number) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string", "", "some_number", "<20", "double")
                .addFuncCallOptionsPattern("some_number", "<20", "<",
                        null, null, 20)
                .build();
        // --input: string some_number = "<=20" {pattern: double}
        FuncCall funcCall13 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string some_number = \"<=20\" {pattern: double}\n"
                        + "MATCH(p:Person) WHERE @some_number(p.some_number) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string", "", "some_number", "<=20", "double")
                .addFuncCallOptionsPattern("some_number", "<=20", "<=",
                        null, null, 20)
                .build();
        // --input: string first_name = 'contains Z' {pattern: string}
        DataFrame expected7 = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"25"}),
                new StringColumn("first_name", new String[]{"Zolly"}),
                new StringColumn("last_name", new String[]{"Wimmers"}),
                new StringColumn("email", new String[]{"zwimmerso@hatena.ne.jp"}),
                new StringColumn("gender", new String[]{"Male"}),
                new StringColumn("ip_address", new String[]{"123.12.225.114/32"}),
                new BoolColumn("bool", new Boolean[]{false}),
                new StringColumn("country", new String[]{"Bosnia and Herzegovina"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "2003-02-12")),
                new FloatColumn("some_number", new Float[]{217.18f}));
        FuncCall funcCall14 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string first_name = 'contains Z' {pattern: string}\n"
                        + "MATCH(p:Person) WHERE @first_name(p.first_name) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string", "", "some_number", ">=975", "double")
                .addFuncParam("string", "","first_name", "contains Z", "string")
                .addFuncCallOptionsPattern("first_name", "contains Z", "contains",
                        null, null, "Z")
                .build();
        // --input: string first_name = 'starts with W' {pattern: string}
        DataFrame expected8 = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"23"}),
                new StringColumn("first_name", new String[]{"Waly"}),
                new StringColumn("last_name", new String[]{"Rogliero"}),
                new StringColumn("email", new String[]{"wroglierom@berkeley.edu"}),
                new StringColumn("gender", new String[]{"Female"}),
                new StringColumn("ip_address", new String[]{"122.90.196.231/32"}),
                new BoolColumn("bool", new Boolean[]{true}),
                new StringColumn("country", new String[]{"Sweden"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "2011-12-18")),
                new FloatColumn("some_number", new Float[]{147.69f}));
        FuncCall funcCall15 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string first_name = 'starts with W' {pattern: string}\n"
                        + "MATCH(p:Person) WHERE @first_name(p.first_name) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string", "","first_name", "starts with W", "string")
                .addFuncCallOptionsPattern("first_name", "starts with W", "starts with",
                        null, null, "W")
                .build();
        // --input: string first_name = 'ends with y' {pattern: string}
        DataFrame expected9 = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"20"}),
                new StringColumn("first_name", new String[]{"Lucius"}),
                new StringColumn("last_name", new String[]{"Edelmann"}),
                new StringColumn("email", new String[]{"ledelmannj@bravesites.com"}),
                new StringColumn("gender", new String[]{"Male"}),
                new StringColumn("ip_address", new String[]{"66.174.30.225/32"}),
                new BoolColumn("bool", new Boolean[]{false}),
                new StringColumn("country", new String[]{"Brazil"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "1999-06-22")),
                new FloatColumn("some_number", new Float[]{378.73f}));
        FuncCall funcCall16 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string first_name = 'ends with s' {pattern: string}\n"
                        + "MATCH(p:Person) WHERE @first_name(p.first_name) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string", "","first_name", "ends with s", "string")
                .addFuncCallOptionsPattern("first_name", "ends with s", "ends with",
                        null, null, "s")
                .build();
        // --input: string country = 'in (Poland, Brazil)' {pattern: string}
        DataFrame expected10 = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"2", "5", "20"}),
                new StringColumn("first_name", new String[]{"Nicholle", "Mitchell", "Lucius",}),
                new StringColumn("last_name", new String[]{"Karoly", "Haglington", "Edelmann"}),
                new StringColumn("email", new String[]{"nkaroly1@alexa.com", "mhaglington4@indiegogo.com",
                        "ledelmannj@bravesites.com"}),
                new StringColumn("gender", new String[]{"Female", "Male", "Male"}),
                new StringColumn("ip_address", new String[]{"255.233.247.118/32", "209.93.181.190/32", "66.174.30.225/32"}),
                new BoolColumn("bool", new Boolean[]{false, true, false}),
                new StringColumn("country", new String[]{"Poland", "Poland", "Brazil"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "2014-02-27",
                                "2020-10-09","1999-06-22")),
                new FloatColumn("some_number", new Float[]{864.09f, 15.22f, 378.73f}));
        FuncCall funcCall17 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string country = 'in (Poland, Brazil)' {pattern: string}\n"
                        + "MATCH(p:Person) WHERE @country(p.country) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string", "", "country", "in (Poland, Brazil)", "string")
                .addFuncCallOptionsPattern("country", "in (Poland, Brazil)", "in",
                        null, null, "Poland", "Brazil")
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
                Arguments.of(Named.of("type: string; operator: in; pattern: string", funcCall17), expected10)
        );
    }

    public static Stream<Arguments> checkRegexSupport_ok() {
        // --input: string email = 'regex ^([A-Za-z0-9_]+@google.com.au)$' {pattern: string}
        DataFrame expected = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"9"}),
                new StringColumn("first_name", new String[]{"Marlie"}),
                new StringColumn("last_name", new String[]{"Mayze"}),
                new StringColumn("email", new String[]{"mmayze8@google.com.au"}),
                new StringColumn("gender", new String[]{"Female"}),
                new StringColumn("ip_address", new String[]{"68.41.25.65/32"}),
                new BoolColumn("bool", new Boolean[]{false}),
                new StringColumn("country", new String[]{"France"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles("yyyy-MM-dd", "2011-11-10")),
                new FloatColumn("some_number", new Float[]{561.72f}));
        FuncCall funcCall = FuncCallBuilder.getBuilder()
                .addQuery("// input: string email = 'regex ^([A-Za-z0-9_]+@google.com.au)$' {pattern: string}\n"
                        + "MATCH(p:Person) WHERE @email(p.email) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string", "", "email", "regex ^([A-Za-z0-9_]+@google.com.au)$", "string")
                .addFuncCallOptionsPattern("email", "regex ^([A-Za-z0-9_]+@google.com.au)$",
                        "regex", null, null, "^([A-Za-z0-9_]+@google.com.au)$")
                .build();
        return Stream.of(Arguments.of(Named.of("type: string; operator: regex; pattern: string",
                funcCall), expected));
    }

    public static Stream<Arguments> checkMultipleParametersSupport_ok() {
        String datePattern = "yyyy-MM-dd";
        // --input: string first_name = "starts with p" {pattern: string}
        //--input: string id = ">1" {pattern :int}
        //--input: bool bool = false
        //--input: string email = "contains com" {pattern: string}
        //--input: string some_number = ">20" {pattern: double}
        //--input: string country = "in (Indonesia)" {pattern: string}
        //--input: string date = "before 1/1/2022" {pattern: datetime}
        DataFrame expected1 = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"13"}),
                new StringColumn("first_name", new String[]{"Pail"}),
                new StringColumn("last_name", new String[]{"Boxell"}),
                new StringColumn("email", new String[]{"pboxellc@moonfruit.com"}),
                new StringColumn("gender", new String[]{"Genderqueer"}),
                new StringColumn("ip_address", new String[]{"2.37.160.155/32"}),
                new BoolColumn("bool", new Boolean[]{false}),
                new StringColumn("country", new String[]{"Indonesia"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "2012-01-14")),
                new FloatColumn("some_number", new Float[]{73.47f}));
        FuncCall funcCall1 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string first_name = \"starts with p\" {pattern: string}\n"
                        + "// input: string id = \">1\" {pattern :int}\n"
                        + "// input: bool bool = false\n"
                        + "// input: string email = \"contains com\" {pattern: string}\n"
                        + "// input: string some_number = \">20\" {pattern: double}\n"
                        + "// input: string country = \"in (Indonesia)\" {pattern: string}\n"
                        + "// input: string date = \"before 1/1/2022\" {pattern: datetime}\n"
                        + "MATCH(p:Person) WHERE @first_name(p.first_name) AND @id(p.id) AND p.bool = @bool AND "
                        + "@email(p.email) AND @some_number(p.some_number) AND @country(p.country) AND @date(p.date)  "
                        + "RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string", "","first_name", "starts with p", "string")
                .addFuncParam("string", "","id", ">1", "int")
                .addFuncParam("bool", "","bool", false, "")
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
                .addQuery("// input: list<string> values\n"
                        + "MATCH(p:Person) WHERE ANY(p.country in @values) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("list", "string", "values", values, "")
                .addFuncCallOptionsPattern("country", "", "",
                        null, null, "Poland", "Brazil")
                .build();
        FuncCall funcCall2 = FuncCallBuilder.getBuilder()
                .addQuery("// input: list<string> values = ['Poland','Brazil']\n"
                        + "MATCH(p:Person) WHERE ANY(p.country in @values) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
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
                new DateTimeColumn("date", new Double[]{parser.parseDateToDouble(datePattern, now.toString())}));
        FuncCall funcCall1 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string date = \"today\" {pattern: datetime}\n"
                        + "MATCH(d:TestDate) WHERE @date(d.date) RETURN d.date AS date;")
                .addFuncParam("string", "","date", "today", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        now.toString(), now.plusDays(1).toString())
                .build();
        // --input: string date = "this week" {pattern: datetime}
        DataFrame expected2 = DataFrame.fromColumns(
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern,
                                now.toString(),
                                dayOfWeek == 1 ? null : yesterday.toString(),
                                lastDayOfWeek.equals(now) ? null : lastDayOfWeek.toString())));
        FuncCall funcCall2 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string date = \"this week\" {pattern: datetime}\n"
                        + "MATCH(d:TestDate) WHERE @date(d.date) RETURN d.date AS date;")
                .addFuncParam("string", "","date", "this week", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        firstDayOfWeek.toString(),
                        lastDayOfWeek.plusDays(1).toString())
                .build();
        // --input: string date = "this month" {pattern: datetime}
        DataFrame expected3 = DataFrame.fromColumns(
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern,
                                now.toString(),
                                dayOfMonth == 1 ? null : yesterday.toString(),
                                lastDayOfWeek.getMonthValue() >  lastDayOfMonth.getMonthValue() || lastDayOfWeek.equals(now)?
                                        null : lastDayOfWeek.toString())));
        FuncCall funcCall3 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string date = \"this month\" {pattern: datetime}\n"
                        + "MATCH(d:TestDate) WHERE @date(d.date) RETURN d.date AS date;")
                .addFuncParam("string", "","date", "this month", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        firstDayOfMonth.toString(),
                        lastDayOfMonth.plusDays(1).toString())
                .build();
        // --input: string date = "this year" {pattern: datetime}
        DataFrame expected4 = DataFrame.fromColumns(
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern,
                                now.toString(),
                                dayOfYear == 1 ? null : yesterday.toString(),
                                lastDayOfWeek.getYear() >  now.getYear() || lastDayOfWeek.equals(now)?
                                        null : lastDayOfWeek.toString(),
                                dayOfLastYear.getYear() == now.getYear() ? dayOfLastYear.toString() : null)));
        FuncCall funcCall4 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string date = \"this year\" {pattern: datetime}\n"
                        + "MATCH(d:TestDate) WHERE @date(d.date) RETURN d.date AS date;")
                .addFuncParam("string","", "date", "this year", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        firstDayOfYear.toString(),
                        lastDayOfYear.plusDays(1).toString())
                .build();
        // --input: string date = "yesterday" {pattern: datetime}
        DataFrame expected5 = DataFrame.fromColumns(
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, yesterday.toString())));
        FuncCall funcCall5 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string date = \"yesterday\" {pattern: datetime}\n"
                        + "MATCH(d:TestDate) WHERE @date(d.date) RETURN d.date AS date;")
                .addFuncParam("string","", "date", "yesterday", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        yesterday.toString(),
                        now.toString())
                .build();
        // --input: string date = "last year" {pattern: datetime}
        DataFrame expected6 = DataFrame.fromColumns(
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern,
                                yesterday.getYear() < now.getYear() ? yesterday.toString() : null,
                                dayOfLastYear.getYear() < now.getYear() ? dayOfLastYear.toString() : null)));
        FuncCall funcCall6 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string date = \"last year\" {pattern: datetime}\n"
                        + "MATCH(d:TestDate) WHERE @date(d.date) RETURN d.date AS date;")
                .addFuncParam("string","","date", "last year", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        firstDayOfYear.minusYears(1).toString(), firstDayOfYear.toString())
                .build();
        // --input: string date = "anytime" {pattern: datetime}
        DataFrame expected7 = DataFrame.fromColumns(
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, now.toString(),
                                yesterday.toString(), lastDayOfWeek.getYear() >  now.getYear() || lastDayOfWeek.equals(now)?
                                        null : lastDayOfWeek.toString(), dayOfLastYear.toString(),
                                "2021-04-09")));
        FuncCall funcCall7 = FuncCallBuilder.getBuilder()
                .addQuery( "// input: string date = \"anytime\" {pattern: datetime}\n"
                        + "MATCH(d:TestDate) WHERE @date(d.date) RETURN d.date AS date;")
                .addFuncParam("string", "","date", "anytime", "datetime")
                .addFuncCallOptionsPattern("date", "", "none", true, true)
                .build();
        // --input: string date = "2021-2022" {pattern: datetime}
        DataFrame expected8 = DataFrame.fromColumns(
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "2021-04-09")));

        FuncCall funcCall8 = FuncCallBuilder.getBuilder()
                .addQuery( "// input: string date = \"2021-2021\" {pattern: datetime}\n"
                        + "MATCH(d:TestDate) WHERE @date(d.date) RETURN d.date AS date;")
                .addFuncParam("string", "","date", "2021-2022", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, true,
                        Year.of(2021).atDay(1).toString(),
                        Year.of(2022).atDay(1).toString())
                .build();
        // --input: string date = "before 1/1/2022" {pattern: datetime}

        FuncCall funcCall9 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string date = \"before 1/1/2022\" {pattern: datetime}\n"
                        + "MATCH(d:TestDate) WHERE @date(d.date) RETURN d.date AS date;")
                .addFuncParam("string", "","date", "before 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("date", "", "before", true, true,
                        Year.of(2022).atDay(1).toString())
                .build();
        // --input: string date = "after 1/1/2022" {pattern: datetime}
        DataFrame expected9 = DataFrame.fromColumns(
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, now.toString(),
                                yesterday.toString(), lastDayOfWeek.toString(), dayOfLastYear.toString())));
        FuncCall funcCall10 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string date = \"after 1/1/2022\" {pattern: datetime}\n"
                        + "MATCH(d:TestDate) WHERE @date(d.date) RETURN d.date AS date;")
                .addFuncParam("string", "","date", "after 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("date", "", "after", true, true,
                        LocalDate.parse("2022-01-01").toString())
                .build();
        // --input: string date = "April 2021" {pattern: datetime}
        FuncCall funcCall11 = FuncCallBuilder.getBuilder()
                .addQuery("// input: string date = \"April 2021\" {pattern: datetime}\n"
                        + "MATCH(d:TestDate) WHERE @date(d.date) RETURN d.date AS date;")
                .addFuncParam("string", "","date", "April 2021", "datetime")
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
//                Arguments.of(Named.of("type: string; operator: last year; pattern: datetime", funcCall6), expected6), // Commented because of NullPointerExc when empty resultSet in Neo4j
                Arguments.of(Named.of("type: string; operator: anytime; pattern: datetime", funcCall7), expected7),
                Arguments.of(Named.of("type: string; operator: range -; pattern: datetime", funcCall8), expected8),
                Arguments.of(Named.of("type: string; operator: before; pattern: datetime", funcCall9), expected8),
                Arguments.of(Named.of("type: string; operator: after; pattern: datetime", funcCall10), expected9),
                Arguments.of(Named.of("type: string; operator: April 2021; pattern: datetime", funcCall11), expected8)
        );
    }

    public static Stream<Arguments> checkAllTypesSupport_ok() {
        DataFrame expected1 = DataFrame.fromColumns(
                new DateTimeColumn("datetime", new Double[]{
                        parser.parseDateToDouble("yyyy-MM-dd'T'HH:mm:ss.SSSX", "2022-06-14T10:02:28.192Z"),
                        parser.parseDateToDouble("yyyy-MM-dd'T'HH:mm:ss.SSSX", "1969-01-12T10:02:28.192Z"),
                        parser.parseDateToDouble("yyyy-MM-dd'T'HH:mm:ss.SSSX", "1909-06-12T10:02:28.192Z")} ),
                new DateTimeColumn("localdatetime", new Double[]{
                        parser.parseDateToDouble("yyyy-MM-dd'T'HH:mm:ss.SSS", "2022-06-14T10:02:30.447"),
                        parser.parseDateToDouble("yyyy-MM-dd'T'HH:mm:ss.SSS", "2001-02-11T10:02:10.445"),
                        parser.parseDateToDouble("yyyy-MM-dd'T'HH:mm:ss.SSS", "1909-01-10T00:02:00.405"),
                }),
                new DateTimeColumn("localtime", new Double[]{
                        parser.parseDateToDouble("HH:mm:ss.SSS", "10:02:31.596"),
                        parser.parseDateToDouble("HH:mm:ss.SSS", "12:02:31.596"),
                        parser.parseDateToDouble("HH:mm:ss.SSS", "00:00:00.000"),
                }),
                new DateTimeColumn("time", new Double[]{
                        parser.parseDateToDouble("yyyy-MM-dd'T'HH:mm:ss.SSS", "1970-01-01T10:02:32.192"),
                        parser.parseDateToDouble("yyyy-MM-dd'T'HH:mm:ss.SSS", "1970-01-01T05:01:59.102"),
                        parser.parseDateToDouble("yyyy-MM-dd'T'HH:mm:ss.SSS", "1970-01-01T00:00:00.000"),
                }),
                new FloatColumn("float", new Float[]{Float.POSITIVE_INFINITY, Float.NEGATIVE_INFINITY, 0.0f}),
                new BigIntColumn("integer", new String[]{"99999999999999", "-132442423112", "-3453463646"}),
                new StringColumn("point", new String[]{"{crs=wgs-84, latitude=13.43, x=56.21, y=13.43, srid=4326, longitude=56.21}",
                "{crs=cartesian-3d, x=0.0, y=4.0, z=1.0, srid=9157}", "{crs=cartesian, x=15.0, y=6.0, srid=7203}"}));
        return Stream.of(Arguments.of(Named.of("ALL TYPES SUPPORT",
                        FuncCallBuilder.fromQuery("MATCH(t:PropertyTest) RETURN t.datetime as datetime, "
                                + "t.localdatetime as localdatetime, t.localtime as localtime, t.time as time, t.float as float, "
                                + "t.integer as integer, t.point as point")), expected1));
    }

    public static Stream<Arguments> checkSupportOfMapReturnType_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new BigIntColumn("p._id", new String[]{"0"}),
                new StringColumn("p._labels", new String[]{"[Person]"}),
                new BigIntColumn("p.id", new String[]{"1"}),
                new StringColumn("p.first_name", new String[]{"Burk"}),
                new StringColumn("p.last_name", new String[]{"Kemery"}),
                new StringColumn("p.email", new String[]{"bkemery0@businesswire.com"}),
                new StringColumn("p.gender", new String[]{"Male"}),
                new StringColumn("p.ip_address", new String[]{"249.64.22.121/32"}),
                new BoolColumn("p.bool", new Boolean[]{true}),
                new StringColumn("p.country", new String[]{"China"}),
                new DateTimeColumn("p.date", parser.parseDatesToDoubles("yyyy-MM-dd", "2017-09-20")),
                new FloatColumn("p.some_number", new Float[]{510.32f}));
        return Stream.of(Arguments.of(Named.of("Map return type support",
                FuncCallBuilder.fromQuery("MATCH(p:Person) RETURN p LIMIT 1")), expected));
    }
}

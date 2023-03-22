package grok_connect.providers.arguments_provider;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.DataFrameBuilder;
import grok_connect.providers.utils.DateParser;
import grok_connect.providers.utils.FuncCallBuilder;
import grok_connect.providers.utils.Parser;
import org.junit.jupiter.api.Named;
import org.junit.jupiter.params.provider.Arguments;
import serialization.*;

import java.util.stream.Stream;

public class Neo4jObjectsMother {
    public static Parser parser = new DateParser();

    public static Stream<Arguments> checkParameterSupport_ok() {
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
                .addQuery("// input: int id = 20\n"
                        + "MATCH(p:Person) WHERE p.id = @id RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("int", "", "id", 20, "")
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
                .addQuery("// input: string first_name = 'starts with W' {pattern: string}\n"
                        + "MATCH(p:Person) WHERE @first_name(p.first_name) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string", "","first_name", "starts with W", "string")
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
                .addQuery("// input: string first_name = 'ends with s' {pattern: string}\n"
                        + "MATCH(p:Person) WHERE @first_name(p.first_name) RETURN p.id AS id, p.first_name AS first_name, "
                        + "p.last_name AS last_name, p.email AS email, p.gender AS gender, p.ip_address AS ip_address, "
                        + "p.bool AS bool, p.country AS country, p.date AS date, p.some_number AS some_number;")
                .addFuncParam("string", "","first_name", "ends with s", "string")
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
        DataFrame expected = DataFrameBuilder.getBuilder()
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
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles("yyyy-MM-dd", "2011-11-10")),
                        "date")
                .setColumn(new FloatColumn(new Float[]{561.72f}), "some_number")
                .build();
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
}

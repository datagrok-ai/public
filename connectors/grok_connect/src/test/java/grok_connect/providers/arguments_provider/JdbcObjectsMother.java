package grok_connect.providers.arguments_provider;

import org.junit.jupiter.api.Named;
import org.junit.jupiter.params.provider.Arguments;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

@SuppressWarnings("unused")
public class JdbcObjectsMother {
    public static Stream<Arguments> test_getParameterNames_ok() {
        String query1 = "--name: getSmth\n" +
                "--connection: SomeConnection\n" +
                "--input: string id = \"1241432\"\n" +
                "SELECT\n" +
                "(select encode(json_data, 'escape')\n" +
                "from some_json_data\n" +
                "where identificator = @id) as \"json\"\n" +
                "(select encode(pdb_data, 'escape')\n" +
                "from pdf_files\n" +
                "where identificator = @id) as \"pdf\",\n" +
                "(select json\n" +
                "from some_sort_of_data\n" +
                "where identificator = @id) as \"nums\",\n" +
                "(select json\n" +
                "from some_json_data\n" +
                "where identificator = @id\n" +
                "limit 1) as \"mocks\";";
        List<String> expectedNames1 = new ArrayList<>();
        expectedNames1.add("id");
        expectedNames1.add("id");
        expectedNames1.add("id");
        expectedNames1.add("id");
        StringBuilder expectedBuffer1 = new StringBuilder("SELECT\n" +
                "(select encode(json_data, 'escape')\n" +
                "from some_json_data\n" +
                "where identificator = ?) as \"json\"\n" +
                "(select encode(pdb_data, 'escape')\n" +
                "from pdf_files\n" +
                "where identificator = ?) as \"pdf\",\n" +
                "(select json\n" +
                "from some_sort_of_data\n" +
                "where identificator = ?) as \"nums\",\n" +
                "(select json\n" +
                "from some_json_data\n" +
                "where identificator = ?\n" +
                "limit 1) as \"mocks\";");
        String query2 = "--input: string first_name = \"starts with p\" {pattern: string}\n"
                + "--input: string id = \">1\" {pattern :int}\n"
                + "--input: bool bool = false\n"
                + "--input: string email = \"contains com\" {pattern: string}\n"
                + "--input: string some_number = \">20\" {pattern: double}\n"
                + "--input: string country = \"in (Indonesia)\" {pattern: string}\n"
                + "--input: string date = \"before 1/1/2022\" {pattern: datetime}\n"
                + "SELECT * FROM mock_data WHERE first_name = @first_name AND id = @id AND bool = @bool "
                + "AND email = @email AND some_number = @some_number "
                + "AND country = @country AND (date >= @date1 AND date <= @date2)\n"
                + "--end";
        List<String> expectedNames2 = new ArrayList<>();
        expectedNames2.add("first_name");
        expectedNames2.add("id");
        expectedNames2.add("bool");
        expectedNames2.add("email");
        expectedNames2.add("some_number");
        expectedNames2.add("country");
        expectedNames2.add("date1");
        expectedNames2.add("date2");
        StringBuilder expectedBuffer2 = new StringBuilder(
                "SELECT * FROM mock_data WHERE first_name = ? AND id = ? AND bool = ? "
                + "AND email = ? AND some_number = ? "
                + "AND country = ? AND (date >= ? AND date <= ?)");
        String query3 = "--input: string first_name = \"starts with p\" {pattern: string}\n"
                + "SELECT * FROM mock_data WHERE first_name = @first_name\n"
                + "--end";
        List<String> expectedNames3 = new ArrayList<>();
        expectedNames3.add("first_name");
        StringBuilder expectedBuffer3 = new StringBuilder("SELECT * FROM mock_data WHERE first_name = ?");
        String query4 =
                "--name: compound activity details for all targets containing @protein\n"
                + "--connection: Chembl\n"
                + "--input: string protein = \"P08172\"\n"
                + "SELECT DISTINCT\n"
                + "  m.chembl_id                      AS compound_chembl_id,\n"
                + "  s.canonical_smiles,\n"
                + "  r.compound_key,\n"
                + "  coalesce(d.pubmed_id::text, d.doi) AS pubmed_id_or_doi,\n"
                + "  a.description                    AS assay_description,\n"
                + "  act.standard_type,\n"
                + "  act.standard_relation,\n"
                + "  act.standard_value,\n"
                + "  act.standard_units,\n"
                + "  act.activity_comment,\n"
                + "  t.chembl_id                      AS target_chembl_id,\n"
                + "  t.pref_name                      AS target_name,\n"
                + "  t.target_type\n"
                + "FROM compound_structures s\n"
                + "  RIGHT JOIN molecule_dictionary m ON s.molregno = m.molregno\n"
                + "  JOIN compound_records r ON m.molregno = r.molregno\n"
                + "  JOIN docs d ON r.doc_id = d.doc_id\n"
                + "  JOIN activities act ON r.record_id = act.record_id\n"
                + "  JOIN assays a ON act.assay_id = a.assay_id\n"
                + "  JOIN target_dictionary t ON a.tid = t.tid\n"
                + "  JOIN target_components tc ON t.tid = tc.tid\n"
                + "  JOIN component_sequences cs ON tc.component_id = cs.component_id\n"
                + "  WHERE cs.accession = @protein;";
        List<String> expectedNames4 = new ArrayList<>();
        expectedNames4.add("protein");
        StringBuilder expectedBuffer4 = new StringBuilder("SELECT DISTINCT\n"
                + "  m.chembl_id                      AS compound_chembl_id,\n"
                + "  s.canonical_smiles,\n"
                + "  r.compound_key,\n"
                + "  coalesce(d.pubmed_id::text, d.doi) AS pubmed_id_or_doi,\n"
                + "  a.description                    AS assay_description,\n"
                + "  act.standard_type,\n"
                + "  act.standard_relation,\n"
                + "  act.standard_value,\n"
                + "  act.standard_units,\n"
                + "  act.activity_comment,\n"
                + "  t.chembl_id                      AS target_chembl_id,\n"
                + "  t.pref_name                      AS target_name,\n"
                + "  t.target_type\n"
                + "FROM compound_structures s\n"
                + "  RIGHT JOIN molecule_dictionary m ON s.molregno = m.molregno\n"
                + "  JOIN compound_records r ON m.molregno = r.molregno\n"
                + "  JOIN docs d ON r.doc_id = d.doc_id\n"
                + "  JOIN activities act ON r.record_id = act.record_id\n"
                + "  JOIN assays a ON act.assay_id = a.assay_id\n"
                + "  JOIN target_dictionary t ON a.tid = t.tid\n"
                + "  JOIN target_components tc ON t.tid = tc.tid\n"
                + "  JOIN component_sequences cs ON tc.component_id = cs.component_id\n"
                + "  WHERE cs.accession = ?;");
        String query5 = "SELECT * FROM MOCK_DATA";
        List<String> expectedNames5 = new ArrayList<>();
        StringBuilder expectedBuffer5 = new StringBuilder(query5);
        String query6 = "SELECT @name FROM MOCK_DATA";
        List<String> expectedNames6 = new ArrayList<>();
        StringBuilder expectedBuffer6 = new StringBuilder(query6);
        String query7 = "-- name: getSomeData\n" +
                "-- friendlyName: getSomeData\n" +
                "-- connection: SomeConnection\n" +
                "-- description: description\n" +
                "-- meta.app: meta\n" +
                "-- meta.model: meta\n" +
                "-- meta.entity: meta\n" +
                "-- meta.ptm: meta\n" +
                "-- input: string ID = 'ID124124'\n" +
                "-- input: string SCHEME = 'some schema'\n" +
                "\n" +
                "-- SET @ID = 'ID124124';\n" +
                "-- set @SCHEMA = 'some schema';\n" +
                "SELECT\n" +
                "    pk.id,\n" +
                "    'Data' as ptm,                            -- comment\n" +
                "    pk.data1,                          -- comment\n" +
                "    pk.data2,                                  -- comment\n" +
                "    pk.data3,                         -- comment\n" +
                "    pk.data4 as data4,                -- comment\n" +
                "    pk.data5,                        -- comment\n" +
                "    pk.data6 as data6,              -- comment\n" +
                "    pk.data7 as data7,              -- comment\n" +
                "    pk.schema as schema             -- comment\n" +
                "FROM mock.mocks as pk\n" +
                "WHERE pk.id = @ID AND pk.schema = @SCHEMA;";
        List<String> expectedNames7 = new ArrayList<>();
        expectedNames7.add("ID");
        expectedNames7.add("SCHEMA");
        StringBuilder expectedBuffer7 = new StringBuilder("SELECT\n" +
                "    pk.id,\n" +
                "    'Data' as ptm,                            -- comment\n" +
                "    pk.data1,                          -- comment\n" +
                "    pk.data2,                                  -- comment\n" +
                "    pk.data3,                         -- comment\n" +
                "    pk.data4 as data4,                -- comment\n" +
                "    pk.data5,                        -- comment\n" +
                "    pk.data6 as data6,              -- comment\n" +
                "    pk.data7 as data7,              -- comment\n" +
                "    pk.schema as schema             -- comment\n" +
                "FROM mock.mocks as pk\n" +
                "WHERE pk.id = ? AND pk.schema = ?;");
        String query8 = "-- name: getData\n" +
                "-- friendlyName: getData\n" +
                "-- connection: SomeConnection\n" +
                "-- description: description\n" +
                "-- tags: tag\n" +
                "-- meta.app: meta\n" +
                "-- meta.entity: entity\n" +
                "-- input: string data1 ='some_text1'\n" +
                "-- input: string data2 ='L'\n" +
                "-- input: int data3 =89\n" +
                "-- input: string data4 ='&'\n" +
                "-- input: int data5 =97\n" +
                "-- input: string data6 ='&'\n" +
                "-- input: string data7 = 'ID124'\n" +
                "\n" +
                "-- SET @data1 = 'some_text1';\n" +
                "-- SET @data2 = 'L';\n" +
                "-- SET @data3 = 89;\n" +
                "-- SET @data4 = '&';\n" +
                "-- SET @data5 = 97;\n" +
                "-- SET @data6 = '&';\n" +
                "-- SET @data7 = 'ID124';\n" +
                "\n" +
                "WITH my_context AS\n" +
                "         (SELECT *,\n" +
                "                 CASE WHEN num2 <> '&' THEN CONCAT(num1, num2) ELSE num1 END as pos\n" +
                "          FROM (SELECT\n" +
                "                    Ar.data7,\n" +
                "                    Ar.data2,\n" +
                "                    Ar.sequence_num,\n" +
                "                    CASE\n" +
                "                        WHEN @data1 = 'word1' THEN Ar.kabat_num1\n" +
                "                        WHEN @data1 = 'word2' THEN Ar.aho_num1\n" +
                "                        WHEN @data1 = 'some_text1' THEN Ar.some_text1_num1\n" +
                "                        WHEN @data1 = 'word3' THEN Ar.imgt_num1\n" +
                "                        WHEN @data1 = 'word4' THEN Ar.martin_num1\n" +
                "                        END as num1,\n" +
                "                    CASE\n" +
                "                        WHEN @data1 = 'word1' THEN Ar.kabat_num2\n" +
                "                        WHEN @data1 = 'word2' THEN Ar.aho_num2\n" +
                "                        WHEN @data1 = 'some_text1' THEN Ar.some_text1_num2\n" +
                "                        WHEN @data1 = 'word3' THEN Ar.imgt_num2\n" +
                "                        WHEN @data1 = 'word4' THEN Ar.martin_num2\n" +
                "                        END as num2,\n" +
                "                    MIN(Ar.aa) as aa\n" +
                "                FROM my_table as Ar\n" +
                "                WHERE Ar.data2 = @data2 AND CASE WHEN @data7 is not NULL THEN Ar.data7 = @data7 ELSE TRUE END\n" +
                "                GROUP BY Ar.data7, Ar.data2, Ar.sequence_num, num1, num2\n" +
                "                ORDER BY num1, num2, Ar.data7) as Ars\n" +
                "          WHERE @data3 <= Ars.num1 AND Ars.num1 <= @data5 AND\n" +
                "                  @data4 <= Ars.num2 AND Ars.num2 <= @data6)\n" +
                "SELECT\n" +
                "    C_data7_pos.data7,\n" +
                "    MIN(Tr2.num1) as start,\n" +
                "    MAX(Tr2.num2) as end,\n" +
                "    GROUP_CONCAT(CASE WHEN Tr2.aa is NULL THEN '-' ELSE Tr2.aa END\n" +
                "                     ORDER BY C_data7_pos.num1, C_data7_pos.some_code SEPARATOR '') as seq\n" +
                "FROM (SELECT\n" +
                "          Cdata7.data7,\n" +
                "          Cpos.num1,\n" +
                "          Cpos.num2,\n" +
                "          Cpos.some_code\n" +
                "      FROM (SELECT\n" +
                "                num1,\n" +
                "                num2,\n" +
                "                some_code\n" +
                "            FROM my_context as Tr\n" +
                "            GROUP BY num1, num2, some_code\n" +
                "            ORDER BY num1, num2) as Cpos\n" +
                "               CROSS JOIN (SELECT\n" +
                "                               Ar.data7\n" +
                "                           FROM my_table as Ar\n" +
                "                           WHERE CASE WHEN @data7 is not NULL THEN Ar.data7 = @data7 ELSE TRUE END\n" +
                "                           GROUP BY Ar.data7) as Cdata7) as C_data7_pos\n" +
                "         LEFT JOIN my_context as Tr2\n" +
                "                   ON Tr2.num1 = C_data7_pos.num1 AND Tr2.num2 = C_data7_pos.num2 AND Tr2.data7 = C_data7_pos.data7\n" +
                "GROUP BY C_data7_pos.data7\n" +
                "UNION\n" +
                "SELECT\n" +
                "    'position_names' as data7,\n" +
                "    -1 as start,\n" +
                "    -1 as end,\n" +
                "    GROUP_CONCAT(schema_pos ORDER BY num1, num2 SEPARATOR ',') as seq\n" +
                "FROM (SELECT\n" +
                "          Tr.num1,\n" +
                "          Tr.num2,\n" +
                "          CONCAT(Tr.num1, ':', Tr.num2) as schema_pos\n" +
                "      FROM my_context as Tr\n" +
                "      GROUP BY num1, num2) as Cpos;";
        List<String> expectedNames8 = new ArrayList<>();
        expectedNames8.add("data1");
        expectedNames8.add("data1");
        expectedNames8.add("data1");
        expectedNames8.add("data1");
        expectedNames8.add("data1");
        expectedNames8.add("data1");
        expectedNames8.add("data1");
        expectedNames8.add("data1");
        expectedNames8.add("data1");
        expectedNames8.add("data1");
        expectedNames8.add("data2");
        expectedNames8.add("data7");
        expectedNames8.add("data7");
        expectedNames8.add("data3");
        expectedNames8.add("data5");
        expectedNames8.add("data4");
        expectedNames8.add("data6");
        expectedNames8.add("data7");
        expectedNames8.add("data7");
        StringBuilder expectedBuffer8 = new StringBuilder("WITH my_context AS\n" +
                "         (SELECT *,\n" +
                "                 CASE WHEN num2 <> '&' THEN CONCAT(num1, num2) ELSE num1 END as pos\n" +
                "          FROM (SELECT\n" +
                "                    Ar.data7,\n" +
                "                    Ar.data2,\n" +
                "                    Ar.sequence_num,\n" +
                "                    CASE\n" +
                "                        WHEN ? = 'word1' THEN Ar.kabat_num1\n" +
                "                        WHEN ? = 'word2' THEN Ar.aho_num1\n" +
                "                        WHEN ? = 'some_text1' THEN Ar.some_text1_num1\n" +
                "                        WHEN ? = 'word3' THEN Ar.imgt_num1\n" +
                "                        WHEN ? = 'word4' THEN Ar.martin_num1\n" +
                "                        END as num1,\n" +
                "                    CASE\n" +
                "                        WHEN ? = 'word1' THEN Ar.kabat_num2\n" +
                "                        WHEN ? = 'word2' THEN Ar.aho_num2\n" +
                "                        WHEN ? = 'some_text1' THEN Ar.some_text1_num2\n" +
                "                        WHEN ? = 'word3' THEN Ar.imgt_num2\n" +
                "                        WHEN ? = 'word4' THEN Ar.martin_num2\n" +
                "                        END as num2,\n" +
                "                    MIN(Ar.aa) as aa\n" +
                "                FROM my_table as Ar\n" +
                "                WHERE Ar.data2 = ? AND CASE WHEN ? is not NULL THEN Ar.data7 = ? ELSE TRUE END\n" +
                "                GROUP BY Ar.data7, Ar.data2, Ar.sequence_num, num1, num2\n" +
                "                ORDER BY num1, num2, Ar.data7) as Ars\n" +
                "          WHERE ? <= Ars.num1 AND Ars.num1 <= ? AND\n" +
                "                  ? <= Ars.num2 AND Ars.num2 <= ?)\n" +
                "SELECT\n" +
                "    C_data7_pos.data7,\n" +
                "    MIN(Tr2.num1) as start,\n" +
                "    MAX(Tr2.num2) as end,\n" +
                "    GROUP_CONCAT(CASE WHEN Tr2.aa is NULL THEN '-' ELSE Tr2.aa END\n" +
                "                     ORDER BY C_data7_pos.num1, C_data7_pos.some_code SEPARATOR '') as seq\n" +
                "FROM (SELECT\n" +
                "          Cdata7.data7,\n" +
                "          Cpos.num1,\n" +
                "          Cpos.num2,\n" +
                "          Cpos.some_code\n" +
                "      FROM (SELECT\n" +
                "                num1,\n" +
                "                num2,\n" +
                "                some_code\n" +
                "            FROM my_context as Tr\n" +
                "            GROUP BY num1, num2, some_code\n" +
                "            ORDER BY num1, num2) as Cpos\n" +
                "               CROSS JOIN (SELECT\n" +
                "                               Ar.data7\n" +
                "                           FROM my_table as Ar\n" +
                "                           WHERE CASE WHEN ? is not NULL THEN Ar.data7 = ? ELSE TRUE END\n" +
                "                           GROUP BY Ar.data7) as Cdata7) as C_data7_pos\n" +
                "         LEFT JOIN my_context as Tr2\n" +
                "                   ON Tr2.num1 = C_data7_pos.num1 AND Tr2.num2 = C_data7_pos.num2 AND Tr2.data7 = C_data7_pos.data7\n" +
                "GROUP BY C_data7_pos.data7\n" +
                "UNION\n" +
                "SELECT\n" +
                "    'position_names' as data7,\n" +
                "    -1 as start,\n" +
                "    -1 as end,\n" +
                "    GROUP_CONCAT(schema_pos ORDER BY num1, num2 SEPARATOR ',') as seq\n" +
                "FROM (SELECT\n" +
                "          Tr.num1,\n" +
                "          Tr.num2,\n" +
                "          CONCAT(Tr.num1, ':', Tr.num2) as schema_pos\n" +
                "      FROM my_context as Tr\n" +
                "      GROUP BY num1, num2) as Cpos;");
        String query9 = "SELECT * FROM MOCK_DATA WHERE name = '--name--';";
        List<String> expectedNames9 = new ArrayList<>();
        StringBuilder expectedBuffer9 = new StringBuilder("SELECT * FROM MOCK_DATA WHERE name = '--name--';");
        String query10 = "-- name: getData\n" +
                "-- friendlyName: getData\n" +
                "-- connection: SomeConnection\n" +
                "-- description: description\n" +
                "-- tags: tag\n" +
                "-- meta.app: meta\n" +
                "-- meta.entity: entity\n"
                +"SELECT * FROM MOCK_DATA WHERE name = '--name--';";
        List<String> expectedNames10 = new ArrayList<>();
        StringBuilder expectedBuffer10 = new StringBuilder("SELECT * FROM MOCK_DATA WHERE name = '--name--';");
        return Stream.of(Arguments.of(Named.of("1 input - 4 parameters", query1),
                expectedNames1, expectedBuffer1),
                Arguments.of(Named.of("6 input - 6 parameters", query2),
                        expectedNames2, expectedBuffer2),
                Arguments.of(Named.of("1 input - 1 parameter", query3),
                        expectedNames3, expectedBuffer3),
                Arguments.of(Named.of("@ in comment, 1 input - 1 parameter", query4),
                        expectedNames4, expectedBuffer4),
                Arguments.of(Named.of("no inputs and parameters", query5),
                        expectedNames5, expectedBuffer5),
                Arguments.of(Named.of("no inputs, but @ in query", query6),
                        expectedNames6, expectedBuffer6),
                Arguments.of(Named.of("2 inputs, 2 parameters, more comments", query7),
                        expectedNames7, expectedBuffer7),
                Arguments.of(Named.of("7 inputs, 7 parameters, more comments", query8),
                        expectedNames8, expectedBuffer8),
                Arguments.of(Named.of("no inputs, no parameters, string as comment", query9),
                        expectedNames9, expectedBuffer9),
                Arguments.of(Named.of("no inputs, no parameters, comments, parameter as a comment", query10),
                        expectedNames10, expectedBuffer10));
    }
}

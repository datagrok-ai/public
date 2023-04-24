package grok_connect.providers.arguments_provider;

import org.junit.jupiter.api.Named;
import org.junit.jupiter.params.provider.Arguments;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

public class JdbcObjectsMother {
    public static Stream<Arguments> test_getParameterNames_ok() {
        String query1 = "--name: get3D\n"
                + "--connection: Mlbdata:MLB\n"
                + "--input: string vid = \"VR000030945\"\n"
                + "SELECT\n"
                + "(select encode(json_data, 'escape')\n"
                + "from db_v2.json_files\n"
                + "where v_id = @vid) as \"json\",\n"
                + "(select encode(pdb_data, 'escape')\n"
                + "from db_v2.pdb_files\n"
                + "where v_id = @vid) as \"pdb\",\n"
                + "(select json\n"
                + "from mlb.jsons_new\n"
                + "where v_id = @vid) as \"real_nums\",\n"
                + "(select json\n"
                + "from db_v2.json_files_observed\n"
                + "where v_id = @vid\n"
                + "limit 1) as \"obs_ptm\";";
        List<String> expectedNames1 = new ArrayList<>();
        expectedNames1.add("vid");
        expectedNames1.add("vid");
        expectedNames1.add("vid");
        expectedNames1.add("vid");
        StringBuilder expectedBuffer1 = new StringBuilder("SELECT\n"
                + "(select encode(json_data, 'escape')\n"
                + "from db_v2.json_files\n"
                + "where v_id = ?) as \"json\",\n"
                + "(select encode(pdb_data, 'escape')\n"
                + "from db_v2.pdb_files\n"
                + "where v_id = ?) as \"pdb\",\n"
                + "(select json\n"
                + "from mlb.jsons_new\n"
                + "where v_id = ?) as \"real_nums\",\n"
                + "(select json\n"
                + "from db_v2.json_files_observed\n"
                + "where v_id = ?\n"
                + "limit 1) as \"obs_ptm\";");
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
                + "AND country = ? AND (date >= ? AND date <= ?)\n"
                + "--end");
        String query3 = "--input: string first_name = \"starts with p\" {pattern: string}\n"
                + "SELECT * FROM mock_data WHERE first_name = @first_name\n"
                + "--end";
        List<String> expectedNames3 = new ArrayList<>();
        expectedNames3.add("first_name");
        StringBuilder expectedBuffer3 = new StringBuilder("SELECT * FROM mock_data WHERE first_name = ?\n"
                        + "--end");
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
                        expectedNames6, expectedBuffer6));
    }
}

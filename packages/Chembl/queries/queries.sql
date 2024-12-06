--name: _protein classification
--friendlyName: Misc | Protein Classification
--connection: Chembl
select protein_class_id, parent_id, pref_name, definition, class_level from protein_classification
--end

--name: _compounds which are selective to one target over a second target
--friendlyName: Browse | Compounds selective to one target over a second target
--connection: Chembl
--input: string selectiveFor = "CHEMBL301"
--input: string over = "CHEMBL4036"
SELECT md.chembl_id,
cs.canonical_smiles
FROM target_dictionary td
  JOIN assays a ON td.tid = a.tid
  JOIN activities act ON a.assay_id = act.assay_id
  JOIN molecule_dictionary md ON md.molregno = act.molregno
  JOIN compound_structures cs ON md.molregno = cs.molregno
    AND act.standard_relation = '='
    AND act.standard_type     IN ('IC50')
    AND act.standard_units    = 'nM'
    AND act.standard_value    < 50
    AND td.chembl_id          = @selectiveFor
INTERSECT
SELECT md.chembl_id,
cs.canonical_smiles
FROM target_dictionary td
  JOIN assays a ON td.tid = a.tid
  JOIN activities act ON a.assay_id = act.assay_id
  JOIN molecule_dictionary md ON md.molregno = act.molregno
  JOIN compound_structures cs ON md.molregno = cs.molregno
AND act.standard_relation     = '='
AND act.standard_type         IN ('IC50')
AND act.standard_units        = 'nM'
AND act.standard_value        > 200
AND td.chembl_id              = @over;
--end

--name: compound activity details for all targets containing @protein
--friendlyName: Browse | Compound activity details for all targets containing @protein
--connection: Chembl
--input: string protein = "P08172"
SELECT DISTINCT
  m.chembl_id                      AS compound_chembl_id,
  s.canonical_smiles,
  r.compound_key,
  coalesce(d.pubmed_id::text, d.doi) AS pubmed_id_or_doi,
  a.description                    AS assay_description,
  act.standard_type,
  act.standard_relation,
  act.standard_value,
  act.standard_units,
  act.activity_comment,
  t.chembl_id                      AS target_chembl_id,
  t.pref_name                      AS target_name,
  t.target_type
FROM compound_structures s
  RIGHT JOIN molecule_dictionary m ON s.molregno = m.molregno
  JOIN compound_records r ON m.molregno = r.molregno
  JOIN docs d ON r.doc_id = d.doc_id
  JOIN activities act ON r.record_id = act.record_id
  JOIN assays a ON act.assay_id = a.assay_id
  JOIN target_dictionary t ON a.tid = t.tid
  JOIN target_components tc ON t.tid = tc.tid
  JOIN component_sequences cs ON tc.component_id = cs.component_id
    AND cs.accession = @protein;
--end

--name: unichemUnitTestQuery
--friendlyName: Misc | Unichem Test
--connection: Unichem
--tags: unit-test
--meta.testExpected: 1
select count(from_id) from src10src11
--end


--name: FracClassification
--friendlyName: Search | FRAC classification
--connection: Chembl 
--input: string level1 = 'MITOSIS AND CELL DIVISION' {choices: Query("SELECT DISTINCT level1_description FROM frac_classification")}
--input: string level2 {nullable: true; choices: Query("SELECT DISTINCT level2_description FROM frac_classification where level1_description = @level1")}
--input: string level3 {nullable: true; choices: Query("SELECT DISTINCT level3_description FROM frac_classification where level2_description = @level2")}
--input: string level4 {nullable: true; choices: Query("SELECT DISTINCT level4_description FROM frac_classification where level3_description = @level3")}
SELECT s.*
FROM compound_structures s
JOIN molecule_frac_classification m
ON s.molregno = m.molregno
JOIN frac_classification f
ON m.frac_class_id = f.frac_class_id
WHERE
  (@level1 is null or @level1 = '' or f.level1_description = @level1) and
  (@level2 is null or @level2 = '' or f.level2_description = @level2) and
  (@level3 is null or @level3 = '' or f.level3_description = @level3) and
  (@level4 is null or @level4 = '' or f.level4_description = @level4)
--end


--name: QueryBySubstructure
--friendlyName: Search | By substructure, country and action type
--connection: Chembl 
--meta.batchMode: true
--input: string substructure = 'c1ccccc1' {semType: Molecule}
--input: string threshold = '0.1' 
--input: string actionType = 'BLOCKER' {choices: Query("SELECT DISTINCT action_type from drug_mechanism")}
--input: string mechanismOfAction = 'Amiloride-sensitive sodium channel, ENaC blocker' {choices: Query("SELECT DISTINCT mechanism_of_action from drug_mechanism where action_type = @actionType")}
--input: string country = 'UK' {choices: Query("SELECT DISTINCT country from research_companies")}
--input: list<string> company = ['GlaxoSmithKline'] {choices: Query("SELECT DISTINCT company from research_companies where country = @country")}
SELECT set_config('rdkit.tanimoto_threshold', @threshold, true);
--batch
SELECT s.*
FROM compound_structures s
INNER JOIN drug_mechanism d
ON s.molregno = d.molregno
INNER JOIN molecule_synonyms m
ON s.molregno = m.molregno
INNER JOIN research_companies r
ON m.res_stem_id = r.res_stem_id
WHERE s.molregno IN (SELECT molregno FROM get_mfp2_neighbors(@substructure))
AND d.action_type = @actionType
AND d.mechanism_of_action = @mechanismOfAction
AND r.country = @country
AND r.company IN (
  SELECT unnest(@company)
)
--end


--name: ByChemblIds
--friendlyName: Search | By ChEMBL ids
--connection: Chembl 
--input: list<string> chemblIds = ['CHEMBL1185', 'CHEMBL1186'] {inputType: TextArea}
SELECT *
FROM molecule_dictionary
WHERE chembl_id IN (
  SELECT unnest(@chemblIds)
);
--end

--name: MolregnoInfo
--connection: Chembl
--tags: panel, widget
--input: string molregno {semType: molregno}
SELECT DISTINCT s.canonical_smiles as smiles, COALESCE(r.country, 'Not found') as country
FROM compound_structures s
LEFT JOIN drug_mechanism d
ON s.molregno = d.molregno
LEFT JOIN molecule_synonyms m
ON s.molregno = m.molregno
LEFT JOIN research_companies r
ON m.res_stem_id = r.res_stem_id
WHERE s.molregno = CAST(@molregno as INTEGER)
--end

--name: ChemblInfo
--connection: Chembl
--tags: panel, widget
--input: string chemblId {semType: CHEMBL_ID}
SELECT DISTINCT s.canonical_smiles as smiles, COALESCE(r.country, 'Not found') as country
FROM molecule_dictionary md
LEFT JOIN compound_structures s
ON md.molregno = s.molregno
LEFT JOIN drug_mechanism d
ON s.molregno = d.molregno
LEFT JOIN molecule_synonyms m
ON s.molregno = m.molregno
LEFT JOIN research_companies r
ON m.res_stem_id = r.res_stem_id
WHERE md.chembl_id = @chemblId
--end
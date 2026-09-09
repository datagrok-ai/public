--name: _protein classification
--friendlyName: Misc | Protein Classification
--description: Returns the complete protein classification hierarchy with preferred names and definitions.
--connection: Chembl
select protein_class_id, parent_id, pref_name, definition, class_level from protein_classification
--end

--name: _compounds which are selective to one target over a second target
--friendlyName: Browse | Compounds Selective For One Target Over Another
--description: Identifies compounds with selective activity, showing high potency for one target and low potency for another target.
--connection: Chembl
--input: string selectiveFor = "CHEMBL301" [ChEMBL target the compound should be potent against]
--input: string over = "CHEMBL4036" [ChEMBL target the compound should be inactive against]
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
--friendlyName: Browse | Compound Activity For Targets Containing Protein
--description: Retrieves compound bioactivity data for all targets containing a specified protein accession identifier.
--connection: Chembl
--input: string protein = "P08172" [UniProt protein accession, e.g. P08172]
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
--description: Test query for Unichem database connectivity and data availability.
--connection: Unichem
--tags: unit-test
--meta.testExpected: 1
select count(from_id) from src10src11
--end


--name: FracClassification
--friendlyName: Search | By FRAC Classification
--description: Searches compound structures by FRAC (Fungicide Resistance Action Committee) mechanism of action.
--connection: Chembl
--input: string mechanism = "tubulin polymerization" {nullable: true; choices: Query("SELECT DISTINCT mechanism_comment FROM pesticide_classification WHERE ref_type = 'FRAC' ORDER BY 1")}
SELECT s.*, p.compound_name, p.mechanism_comment
FROM compound_structures s
JOIN pesticide_class_mapping m
ON s.molregno = m.molregno
JOIN pesticide_classification p
ON m.pest_class_id = p.pest_class_id
WHERE p.ref_type = 'FRAC'
  AND (@mechanism is null or @mechanism = '' or p.mechanism_comment = @mechanism)
--end


--name: QueryBySubstructure
--friendlyName: Search | By Substructure And Action Type
--description: Complex search combining molecular similarity with drug mechanism action type.
--connection: Chembl
--meta.batchMode: true
--input: string substructure = 'c1ccccc1' {semType: Molecule}
--input: string threshold = '0.1'
--input: string actionType = 'BLOCKER' {choices: Query("SELECT DISTINCT action_type from drug_mechanism")}
--input: string mechanismOfAction = 'Amiloride-sensitive sodium channel, ENaC blocker' {choices: Query("SELECT DISTINCT mechanism_of_action from drug_mechanism where action_type = @actionType")}
SELECT set_config('rdkit.tanimoto_threshold', @threshold, true);
--batch
SELECT s.*
FROM compound_structures s
INNER JOIN drug_mechanism d
ON s.molregno = d.molregno
WHERE s.molregno IN (SELECT molregno FROM get_mfp2_neighbors(@substructure))
AND d.action_type = @actionType
AND d.mechanism_of_action = @mechanismOfAction
--end


--name: ByChemblIds
--friendlyName: Search | By ChEMBL IDs
--description: Retrieves molecule dictionary information for a list of ChEMBL identifiers.
--connection: Chembl
--input: list<string> chemblIds = ['CHEMBL1185', 'CHEMBL1186'] {inputType: TextArea} [List of public ChEMBL compound identifiers]
SELECT *
FROM molecule_dictionary
WHERE chembl_id IN (
  SELECT unnest(@chemblIds)
);
--end

--name: MolregnoInfo
--friendlyName: Misc | Compound Info by Molregno
--description: Retrieves compound SMILES, preferred name and max clinical phase for a given molregno identifier.
--connection: Chembl
--tags: panel, widget
--input: int molregno {semType: molregno} [ChEMBL internal compound registration number]
SELECT s.canonical_smiles as smiles, COALESCE(md.pref_name, 'Not found') as name, md.max_phase
FROM compound_structures s
JOIN molecule_dictionary md
ON s.molregno = md.molregno
WHERE s.molregno = CAST(@molregno as INTEGER)
--end

--name: ChemblInfo
--friendlyName: Misc | Compound Info by ChEMBL ID
--description: Retrieves compound SMILES, preferred name and max clinical phase for a given ChEMBL identifier.
--connection: Chembl
--tags: panel, widget
--input: string chemblId {semType: CHEMBL_ID} [Public ChEMBL compound identifier, e.g. CHEMBL1185]
SELECT s.canonical_smiles as smiles, COALESCE(md.pref_name, 'Not found') as name, md.max_phase
FROM molecule_dictionary md
LEFT JOIN compound_structures s
ON md.molregno = s.molregno
WHERE md.chembl_id = @chemblId
--end


--name: FracClassificationWithSubstructure
--friendlyName: Search | By FRAC Classification And Substructure
--description: Combines FRAC mechanism of action search with molecular substructure matching.
--connection: Chembl
--input: string mechanism = "C14-demethylase in sterol biosynthesis (erg11/cyp51)" {nullable: true; choices: Query("SELECT DISTINCT mechanism_comment FROM pesticide_classification WHERE ref_type = 'FRAC' ORDER BY 1")}
--input: string substructure = "Clc1ccccc1" {semType: Substructure}
SELECT s.*, p.compound_name, p.mechanism_comment
FROM compound_structures s
JOIN pesticide_class_mapping m
ON s.molregno = m.molregno
JOIN pesticide_classification p
ON m.pest_class_id = p.pest_class_id
WHERE p.ref_type = 'FRAC'
  AND (@mechanism is null or @mechanism = '' or p.mechanism_comment = @mechanism)
  AND s.canonical_smiles::mol @> @substructure::qmol
--end
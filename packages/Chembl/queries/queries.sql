--name: protein classification
--connection: Chembl
select protein_class_id, parent_id, pref_name, definition, class_level from protein_classification
--end

--name: bioactivity data for bacterial targets for @organism
--connection: Chembl
--input: string organism = "Shigella" {suggestions: Chembl:organisms}
--tags: unit-test
SELECT md.chembl_id AS compound_chembl_id,
cs.canonical_smiles,
act.standard_type,
act.standard_value,
act.standard_units,
td.chembl_id AS target_chembl_id,
td.organism,   td.pref_name
FROM target_dictionary td
  JOIN assays a ON td.tid = a.tid
  JOIN activities act ON a.assay_id = act.assay_id
  JOIN molecule_dictionary md ON act.molregno = md.molregno
  JOIN compound_structures cs ON md.molregno   = cs.molregno
  JOIN organism_class oc ON td.tax_id = oc.tax_id
    AND td.organism = @organism
    AND oc.L1 = 'Bacteria';
--end

--name: activity details for compound and all its salts which have an IC50 bioactivity value in nM against a target of interest
--connection: Chembl
--input: string compound = "CHEMBL192"
--input: string target = "CHEMBL1827"
SELECT m.chembl_id AS compound_chembl_id,
s.canonical_smiles,
r.compound_key,
coalesce(d.pubmed_id::text,d.doi) AS pubmed_id_or_doi,
a.description                   AS assay_description,
act.standard_type,
act.standard_relation,
act.standard_value,
act.standard_units,
act.activity_comment,
t.chembl_id                    AS target_chembl_id,
t.pref_name                    AS target_name,
t.organism                     AS target_organism
FROM compound_structures s
  RIGHT JOIN molecule_dictionary m ON s.molregno = m.molregno
  JOIN compound_records r ON m.molregno = r.molregno
  JOIN docs d ON r.doc_id = d.doc_id
  JOIN activities act ON r.record_id = act.record_id
  JOIN assays a ON act.assay_id = a.assay_id
  JOIN target_dictionary t ON a.tid = t.tid
    AND t.chembl_id      = @target
    AND m.chembl_id IN
         (SELECT DISTINCT
            m1.chembl_id
          FROM molecule_dictionary m1
            JOIN molecule_hierarchy mh ON mh.molregno = m1.molregno
            JOIN molecule_dictionary m2 ON mh.parent_molregno = m2.molregno
              AND m2.chembl_id = @compound)
    AND act.standard_type = 'IC50'
    AND act.standard_units = 'nM';
--end

--name: compounds which are selective to one target over a second target
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

--name: target ChEMBL_ID, target_name, target_type, protein accessions and sequences for all protein targets
--connection: Chembl
SELECT t.chembl_id AS target_chembl_id,
t.pref_name        AS target_name,
t.target_type,
c.accession        AS protein_accession,
c.sequence         AS protein_sequence
FROM target_dictionary t
  JOIN target_type tt ON t.target_type = tt.target_type
  JOIN target_components tc ON t.tid = tc.tid
  JOIN component_sequences c ON tc.component_id = c.component_id
AND tt.parent_type  = 'PROTEIN';
--end

--name: PK data from 'Curated Drug Pharmacokinetic Data' source for @drug
--connection: Chembl
--input: string drug = "LEVOFLOXACIN"
SELECT DISTINCT
  d.title,
  min(case when ap.standard_type = 'DATASET' then coalesce(ap.standard_value::text, ap.standard_text_value) else null end)         dataset,
  a.assay_id,
  a.description,
  min(case when actp.standard_type = 'DOSED_COMPOUND_NAME' then coalesce(actp.standard_value::text, actp.standard_text_value) || ' ' || actp.standard_units  else null end)  dosed_compound_name,
  min(case when actp.standard_type = 'DOSE'                then coalesce(actp.standard_value::text, actp.standard_text_value) || ' ' || actp.standard_units  else null end)  dose,
  min(case when actp.standard_type = 'DOSAGE_FORM'         then coalesce(actp.standard_value::text, actp.standard_text_value) || ' ' || actp.standard_units  else null end)  dosage_form,
  min(case when actp.standard_type = 'REGIMEN'             then coalesce(actp.standard_value::text, actp.standard_text_value) || ' ' || actp.standard_units  else null end)  regimen,
  min(case when actp.standard_type = 'ROUTE'               then coalesce(actp.standard_value::text, actp.standard_text_value)  else null end) route,
  min(case when actp.standard_type = 'GENDER'              then coalesce(actp.standard_value::text, actp.standard_text_value)  else null end) gender,
  min(case when actp.standard_type = 'AGE_RANGE'           then coalesce(actp.standard_value::text, actp.standard_text_value)  else null end) age_range,
  min(case when actp.standard_type = 'HEALTH_STATUS'       then coalesce(actp.standard_value::text, actp.standard_text_value)  else null end) health_status,
  min(case when actp.standard_type = 'TISSUE'              then coalesce(actp.standard_value::text, actp.standard_text_value)  else null end) tissue,
  cr.molregno,
  cr.compound_name,
  act.activity_id,
  act.toid,
  act.standard_type,
  act.standard_relation,
  act.standard_value,
  act.standard_units,
  act.activity_comment
FROM source s
  JOIN compound_records cr ON s.src_id = cr.src_id
  JOIN docs d ON d.doc_id = cr.doc_id
  JOIN activities act ON cr.record_id = act.record_id AND cr.doc_id = act.doc_id
  JOIN activity_properties actp ON act.activity_id = actp.activity_id
  JOIN assays a ON act.assay_id = a.assay_id
  JOIN assay_parameters ap ON a.assay_id = ap.assay_id
                              AND s.src_description = 'Curated Drug Pharmacokinetic Data'
                              AND cr.compound_name = @drug
GROUP BY d.title, a.assay_id, a.description, cr.molregno, cr.compound_name, act.activity_id, act.toid,
  act.standard_type, act.standard_relation, act.standard_value, act.standard_units, act.activity_comment
ORDER BY cr.compound_name, act.toid, act.standard_type;
--end

--name: compound activity details for all targets containing @protein
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

--name: compound activity details for @target
--connection: Chembl
--input: string target = "CHEMBL1827"
--meta.cache: true
--meta.localCache: true
--meta.invalidate: 0 0 * ? * * *
SELECT m.chembl_id AS compound_chembl_id,
s.canonical_smiles,
r.compound_key,
coalesce(d.pubmed_id::text, d.doi) AS pubmed_id_or_doi,
a.description                   AS assay_description,   act.standard_type,
act.standard_relation,
act.standard_value,
act.standard_units,
act.activity_comment
FROM compound_structures s,
molecule_dictionary m,
compound_records r,
docs d,
activities act,
assays a,
target_dictionary t
WHERE s.molregno     = m.molregno
AND m.molregno       = r.molregno
AND r.record_id      = act.record_id
AND r.doc_id         = d.doc_id
AND act.assay_id     = a.assay_id
AND a.tid            = t.tid
AND t.chembl_id      = @target;
--end

--name: all chembl ids with inchi keys
--connection: Chembl
select
  molecule_dictionary.chembl_id,
  compound_structures.standard_inchi_key
from
  molecule_dictionary
  left join compound_structures on molecule_dictionary.molregno = compound_structures.molregno
--end

--name: allChemblStructures
--connection: Chembl
select
  canonical_smiles
from
  compound_structures
--end

--name: unichemUnitTestQuery
--connection: Unichem
--tags: unit-test
--meta.testExpected: 1
select count(from_id) from src10src11
--end

--name: _compoundNames
--connection: Chembl
--input: string sub
--tags: unit-test
select distinct compound_name from compound_structures
where compound_name ilike '%' || @sub || '%'
limit 50
--end

--name: _organisms
--connection: Chembl
--input: string sub
--tags: unit-test
select distinct organism from target_dictionary
where organism ilike '%' || @sub || '%'
limit 50
--end

--name: _proteinTypes
--connection: Chembl
--input: string sub
--tags: unit-test
select distinct short_name from protein_classification
where short_name ilike '%' || @sub || '%'
limit 50
--end

--name: _targetTypes
--connection: Chembl
--input: string sub
select target_type from target_type
where target_type ilike '%' || @sub || '%'
limit 50
--end

--name: _assayTypes
--connection: Chembl
--input: string sub
select assay_type from assay_type
where assay_type ilike '%' || @sub || '%'
limit 50
--end

--name: _relationshipTypes
--connection: Chembl
--input: string sub
select relationship_type from relationship_type
where relationship_type ilike '%' || @sub || '%'
limit 50
--end

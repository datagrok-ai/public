--name: cbAllChemblStructures
--connection: Chembl
select
  canonical_smiles, molregno
from
  compound_structures
limit 1000
--end

--name: cbFindByMolregno
--connection: Chembl
--input: int molregno
select
  compound_structures.canonical_smiles,
  compound_structures.molregno,
  compound_records.compound_name,
  compound_properties.num_ro5_violations,
  molecule_synonyms.synonyms,
  molecule_dictionary.max_phase,
  molecule_dictionary.chembl_id
from
  compound_structures
  left join compound_records on compound_records.molregno = compound_structures.molregno
  left join molecule_dictionary on molecule_dictionary.molregno = compound_structures.molregno
  left join molecule_synonyms on molecule_synonyms.molregno = compound_structures.molregno
  left join compound_properties on compound_properties.molregno = compound_structures.molregno
where
 compound_records.molregno = @molregno
--end

--name: cbChemblBrowserQuery
--connection: Chembl
--input: string substructure
--input: string molecule_type
--input: string subname
--input: int max_phase
--input: int num_ro5_violations
select
  compound_structures.canonical_smiles,
  compound_structures.molregno,
  compound_records.compound_name,
  compound_properties.num_ro5_violations,
  molecule_synonyms.synonyms,
  molecule_dictionary.molecule_type,
  molecule_dictionary.max_phase,
  molecule_dictionary.chembl_id
from
  compound_structures
  left join compound_records on compound_records.molregno = compound_structures.molregno
  left join molecule_dictionary on molecule_dictionary.molregno = compound_structures.molregno
  left join molecule_synonyms on molecule_synonyms.molregno = compound_structures.molregno
  left join compound_properties on compound_properties.molregno = compound_structures.molregno
where
 (compound_structures.canonical_smiles ilike '%' || @substructure || '%' or @substructure is null)
 and (compound_records.compound_name ilike '%' || @subname || '%' or @subname is null)
 and (molecule_dictionary.max_phase = @max_phase or @max_phase = -1)
 and (compound_properties.num_ro5_violations = @num_ro5_violations or @num_ro5_violations = -1)
 and (molecule_dictionary.molecule_type = @molecule_type or @molecule_type = 'All')
 limit 10000
--end

--name: compoundsSelectiveToOneTargetOverSecond
--connection: Chembl
--input: string selectiveFor = CHEMBL301
--input: string over = CHEMBL4036
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

--name: compoundActivityDetailsForTarget
--connection: Chembl
--input: string target = CHEMBL1827
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

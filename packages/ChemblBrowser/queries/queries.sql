--name: allChemblStructures
--connection: chembl:chembl
select
  canonical_smiles, molregno, standard_inchi
from
  compound_structures
limit 100
--end

--name: proteinClassification
--connection: chembl:chembl
select pref_name, definition, class_level from protein_classification
--end


--name: allChemblIdsWithInchiKeys
--connection: chembl:chembl
select
  molecule_dictionary.chembl_id,
  compound_structures.standard_inchi_key
from
  molecule_dictionary
  left join compound_structures on molecule_dictionary.molregno = compound_structures.molregno
limit 1000
--end

--name: CompoundProperties
--connection: chembl:chembl
--input: int molregno
select * from
  compound_properties
where molregno = @molregno
--end



--name: compoundsSelectiveToOneTargetOverSecond
--connection: chembl:chembl
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
--connection: chembl:chembl
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
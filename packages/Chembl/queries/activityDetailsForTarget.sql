--name: compound activity details for @target
--friendlyName: Browse | Compound activity details for @target
--connection: Chembl
--input: string target = "CHEMBL1827"
--meta.searchPattern: "compound activity details for target ${target}"
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
AND t.chembl_id      ILIKE @target;
--end
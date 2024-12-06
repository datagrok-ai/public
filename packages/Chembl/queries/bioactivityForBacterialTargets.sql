--name: _bioactivity data for bacterial targets for @organism
--friendlyName: Browse | Bioactivity for bacterial targets for @organism
--connection: Chembl
--input: string organism = "Shigella" {suggestions: Chembl:organisms}
--meta.searchPattern: "Bioactivity for bacterial targets for ${organism}"
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
    AND td.organism ILIKE @organism
    AND oc.L1 = 'Bacteria';
--end

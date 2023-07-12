--name: _cbAllChemblStructures
--friendlyName: Browse | All ChEMBL structures
--connection: Chembl
select
  canonical_smiles, molregno
from
  compound_structures
limit 1000
--end

--name: _cbAllChemblNumberOfStructures
--friendlyName: Browse | All ChEMBL structures
--input: int maxNumberOfMolecules = 1000
--connection: Chembl
select
  canonical_smiles, molregno
from
  compound_structures
limit @maxNumberOfMolecules
--end

--name: _cbOfStructuresByOrganism
--friendlyName: Bioactivity for bacterial targets
--input: int maxNumberOfMolecules = 1000
--input: string organism = "Shigella" {suggestions: Chembl:organisms}
--connection: Chembl
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
    AND oc.L1 = 'Bacteria'
limit @maxNumberOfMolecules;
--end

--name: _cbFindByMolregno
--friendlyName: Misc | Find by Molregno
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


--name: _cbChemblBrowserQuery
--friendlyName: Browse | Summary
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
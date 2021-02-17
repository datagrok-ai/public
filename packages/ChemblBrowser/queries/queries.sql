--name: allChemblStructures
--connection: chembl:chembl
select
  canonical_smiles, molregno
from
  compound_structures
limit 40
--end

--name: proteinClassification
--connection: chembl:chembl
select protein_class_id, parent_id, pref_name, definition, class_level from protein_classification
--end


--name: allChemblIdsWithInchiKeys
--connection: chembl:chembl
select
  molecule_dictionary.chembl_id,
  compound_structures.standard_inchi_key
from
  molecule_dictionary
  left join compound_structures on molecule_dictionary.molregno = compound_structures.molregno
--end




--name: allChemblStructures
--connection: chembl:chembl
select
  canonical_smiles, molregno
from
  compound_structures
limit 1000
--end


--name: FindBySubstructure
--connection: chembl:chembl
--input: string sub
select
  canonical_smiles, molregno
from
  compound_structures
where
  canonical_smiles ilike '%' || @sub || '%'
limit 1000
--end

--name: FindByName
--connection: chembl:chembl
--input: string sub
SELECT compound_records.molregno,  compound_structures.canonical_smiles, compound_records.compound_name
FROM compound_records
INNER JOIN compound_structures
ON compound_records.molregno = compound_structures.molregno
WHERE compound_name ilike '%' || @sub || '%'
limit 1000
--end


--name: FindByMolregno
--connection: chembl:chembl
--input: int molregno
select
  canonical_smiles, molregno
from
  compound_structures
where molregno = @molregno
limit 1000
--end

--name: FindByRO5
--connection: chembl:chembl
--input: int num_ro5_violations
SELECT compound_structures.molregno, compound_structures.canonical_smiles, compound_properties.num_ro5_violations
FROM compound_structures
INNER JOIN compound_properties
ON compound_structures.molregno = compound_properties.molregno
WHERE num_ro5_violations = @num_ro5_violations
limit 1000
--end

--name: ShowSynonyms
--connection: chembl:chembl
SELECT compound_structures.molregno, compound_structures.canonical_smiles, molecule_synonyms.synonyms
FROM compound_structures
INNER JOIN molecule_synonyms
ON compound_structures.molregno = molecule_synonyms.molregno
limit 1000
--end

--name: ShowChemblID
--connection: chembl:chembl
SELECT compound_structures.molregno, compound_structures.canonical_smiles, molecule_dictionary.chembl_id
FROM compound_structures
INNER JOIN molecule_dictionary
ON compound_structures.molregno = molecule_dictionary.molregno
limit 1000
--end

--name: FilterByMoleculeType
--connection: chembl:chembl
--input: string molecule_type
SELECT compound_structures.molregno, compound_structures.canonical_smiles, molecule_dictionary.molecule_type
FROM compound_structures
INNER JOIN molecule_dictionary
ON compound_structures.molregno = molecule_dictionary.molregno
WHERE molecule_type = @molecule_type
limit 1000
--end

--name: FilterByMaxPhase
--connection: chembl:chembl
--input: int max_phase
SELECT compound_structures.molregno, compound_structures.canonical_smiles, molecule_dictionary.max_phase
FROM compound_structures
INNER JOIN molecule_dictionary
ON compound_structures.molregno = molecule_dictionary.molregno
WHERE max_phase = @max_phase
limit 1000
--end


--name: ChemblBrowserQuery
--connection: chembl:chembl
--input: string substructure = c1
--input: string subname = benzil
--input: int max_phase = 0
--input: string molecule_type = Protein { choices: Query("SELECT DISTINCT molecule_type FROM molecule_dictionary") }
--input: int num_ro5_violations = 0
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
 and (molecule_dictionary.max_phase = @max_phase or @max_phase is null)
 and (compound_properties.num_ro5_violations = @num_ro5_violations or @num_ro5_violations is null)
limit 1000
--end



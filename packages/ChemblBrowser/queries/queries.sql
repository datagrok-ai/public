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




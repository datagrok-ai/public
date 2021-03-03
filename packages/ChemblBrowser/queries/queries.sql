--name: allChemblStructures
--connection: chembl:chembl
select
  canonical_smiles, molregno
from
  compound_structures
limit 100
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
--end

--name: FindByName
--connection: chembl:chembl
--input: string sub
SELECT *
FROM compound_records
WHERE compound_name ilike '%' || @sub || '%'
LEFT JOIN compound_structures
ON compound_records.molregno = compound_structures.molregno
--end


--name: FindByMolregno
--connection: chembl:chembl
--input: int molregno
select
  canonical_smiles, molregno
from
  compound_structures
where molregno = @molregno
--end



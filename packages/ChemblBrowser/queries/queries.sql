--name: allChemblStructures
--input: int molregno
--connection: chembl:chembl
select
  molregno, canonical_smiles
from
  compound_structures
where molregno=@molregno
--end

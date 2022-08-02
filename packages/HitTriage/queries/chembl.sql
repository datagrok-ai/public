--name: bySmiles
--connection: chembl:chembl
--input: int molregno
select molregno, canonical_smiles from compound_structures where molregno = @molregno
--end


--name: Top100Smiles
--connection: chembl:chembl
select molregno, canonical_smiles from compound_structures limit 100
--end


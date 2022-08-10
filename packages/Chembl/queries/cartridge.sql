--name: patternSimilaritySearch
--connection: Chembl
--input: string pattern {semType: Molecule}
--input: int maxRows = 1000
select fps.molregno, cs.canonical_smiles as smiles, tanimoto_sml(morganbv_fp(mol_from_smiles($1::cstring)),mfp2) as similarity
from rdk.fps fps
join compound_structures cs on cs.molregno = fps.molregno
where morganbv_fp(@pattern::mol)%mfp2
order by morganbv_fp(mol_from_smiles($1::cstring))<%>mfp2
limit @maxRows
--end

--name: patternSimilaritySearchWithThreshold
--connection: Chembl
--meta.batchMode: true
--input: string pattern {semType: Molecule}
--input: string threshold = "0.6"
select set_config('rdkit.tanimoto_threshold', @threshold, true);
--batch
select * from get_mfp2_neighbors(@pattern);
--end

--name: patternSubstructureSearch
--connection: Chembl
--input: string pattern {semType: Substructure}
--input: int maxRows = 1000
 select molregno,m as smiles from rdk.mols where m@>@pattern::qmol
 limit @maxRows
--end
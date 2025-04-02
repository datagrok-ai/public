--name: patternSimilaritySearch
--friendlyName: Search | Pattern Similarity
--connection: Chembl
--input: string pattern {semType: Molecule}
--input: int maxRows = 1000
select md.chembl_id, fps.molregno, cs.canonical_smiles as smiles, tanimoto_sml(morganbv_fp(mol_from_smiles($1::cstring)),mfp2) as similarity
from rdk.fps fps
join compound_structures cs on cs.molregno = fps.molregno
join molecule_dictionary md on cs.molregno = md.molregno
where morganbv_fp(@pattern::mol)%mfp2
order by morganbv_fp(mol_from_smiles($1::cstring))<%>mfp2
limit @maxRows
--end

--name: patternSimilaritySearchWithThreshold
--friendlyName: Search | Pattern Similarity With Threshold
--meta.description: Search for a given pattern in the ChEMBL database with a specified threshold of similarity.
--connection: Chembl
--meta.batchMode: true
--input: string pattern {semType: Molecule}
--input: double threshold = 0.6 { min: 0, max: 1 }
select set_config('rdkit.tanimoto_threshold', @threshold::text, true);
--batch
select molregno, m as molecule, similarity from get_mfp2_neighbors(@pattern);
--end

--name: patternSubstructureSearch
--friendlyName: Search | Substructure
--meta.description: Search for a given substructure in the ChEMBL database.
--connection: Chembl
--input: string pattern {semType: Substructure}
--input: int maxRows = 1000
 select md.chembl_id, md.molregno, m as smiles from rdk.mols mols
 join molecule_dictionary md on mols.molregno = md.molregno
 where m@>@pattern::qmol
 limit @maxRows
--end

--name: ChemblNumberOfStructures
--friendlyName: Browse | Specified number of ChEMBL structures
--input: int maxNumberOfMolecules = 1000
--connection: Chembl
select
  canonical_smiles, molregno
from
  compound_structures
limit @maxNumberOfMolecules
--end

--name: ChemblMolregNoBySmiles
--friendlyName: Chembl Molregno by smiles
--input: string smiles {semType: Molecule}
--connection: Chembl
select
  molregno
from
  compound_structures
where
  canonical_smiles = @smiles
limit 1
--end

--name: StructuresByOrganism
--friendlyName: Chembl Targets by organism
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
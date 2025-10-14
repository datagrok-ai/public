--name: getCompounds
--connection: moltrack
SELECT id, canonical_smiles from moltrack.compounds;
--end
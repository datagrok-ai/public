/* Dataset aliases → platform locations. `{dataset}` also accepts a literal `System:…` path. */
import {dataset} from '../../src/registry.js';

dataset('spgi', {path: 'System:AppData/Chem/tests/spgi-100.csv', aliases: ['spgi-100'],
  description: 'SMILES + numeric activity, 100 rows'});
dataset('mol1K', {path: 'System:AppData/Chem/mol1K.csv', aliases: ['mol1k'],
  description: '1000 molecules with pIC50_HIV_Integrase and Q (comes with the published Chem package)'});
dataset('smiles', {path: 'System:DemoFiles/chem/smiles.csv',
  description: '1000 drug-like molecules in canonical_smiles'});
dataset('FASTA_PT_activity', {path: 'System:AppData/Bio/samples/FASTA_PT_activity.csv', aliases: ['fasta-pt-activity', 'peptides with activity'],
  description: '99 peptides: cluster, sequence_id, sequence (16-mers), activity, is_cliff (comes with the published Bio package)'});
dataset('demog', {path: 'System:DemoFiles/demog.csv', description: 'the demographics demo table'});
dataset('demog-1000', {path: 'System:DemoFiles/demog-1000.csv',
  description: 'a stratified 1000-row subset of demog (same SEX / RACE / DIS_POP proportions) — the table for viewer features: every paint costs one marker per row'});
dataset('cars', {path: 'System:DemoFiles/cars.csv'});
dataset('iris', {path: 'System:DemoFiles/iris.csv', description: '150 flowers: four measurements and the Species category'});
dataset('earthquakes', {path: 'System:DemoFiles/geo/earthquakes.csv',
  description: '2426 quakes with Latitude / Longitude / Depth / Magnitude — the geo table the map viewers bind to'});
dataset('beer', {path: 'System:DemoFiles/beer.csv', description: '118 beers, 33 columns; Aroma is a long-text column, so its default filter is a text filter'});
dataset('curves', {path: 'System:DemoFiles/curves.csv', description: 'fit curves ("multiple prefit" carries the fit semantic type) next to a smiles column'});
dataset('smiles', {path: 'System:DemoFiles/chem/smiles.csv',
  description: '1000 ChEMBL molecules: molregno, canonical_smiles (a Molecule column the Chem package renders and offers its Current Value actions on) and RDKit descriptors'});
dataset('helm-peptides', {path: 'System:DemoFiles/chem/peptides/HELM.csv', aliases: ['helm'],
  description: '540 peptides in HELM notation with Activity — a Macromolecule column the Helm package renders and edits'});
dataset('spgi-linked1', {path: 'System:AppData/ApiTests/datasets/SPGI-linked1.csv',
  description: 'the table linked to spgi-100 by Id / Concept Id (the ApiTests package must be published)'});
dataset('spgi-linked2', {path: 'System:AppData/ApiTests/datasets/SPGI-linked2.csv',
  description: 'the table linked to spgi-linked1 by four key columns (the ApiTests package must be published)'});

/* The datasets Chem's features open: the package's own test files on the stand
   (`System:AppData/Chem/...`, published with the package). */
import {dataset} from '@datagrok-libraries/bdd';

dataset('smiles-50', {path: 'System:AppData/Chem/tests/smiles-50.csv',
  description: '50 molecules in canonical_smiles'});
dataset('mol1K.sdf', {path: 'System:AppData/Chem/mol1K.sdf', aliases: ['mol1K sdf'],
  description: '1000 molecules as V2000 molblocks with their properties'});
dataset('ApprovedDrugs2015', {path: 'System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf', aliases: ['approved drugs sdf'],
  description: 'approved drugs of 2015 as V3000 molblocks'});
dataset('sar-small', {path: 'System:DemoFiles/chem/sar_small.csv', aliases: ['sar_small'],
  description: '200 molecules of one SAR series in "smiles", with activity'});
dataset('molecules.mol2', {path: 'System:AppData/Chem/molecules.mol2', aliases: ['molecules mol2'],
  description: 'three TRIPOS molecule blocks'});
dataset('aspirin.mol', {path: 'System:AppData/Chem/tests/molfileV2000.mol', aliases: ['molfile v2000'],
  description: 'one V2000 molblock'});
dataset('test-reactions', {path: 'System:AppData/Chem/test-reactions.csv', aliases: ['test_reactions'],
  description: 'reaction SMILES in a "reaction" column'});
dataset('test_mixtures', {path: 'System:AppData/Chem/test_mixtures.csv', aliases: ['test-mixtures'],
  description: 'chemical mixtures as Mixfile JSON in a "mixture" column'});
dataset('mmp-demo', {path: 'System:DemoFiles/chem/mmp_demo.csv', aliases: ['mmp_demo'],
  description: 'molecules with CYP3A4 and hERG_pIC50 activities'});
dataset('smiles-2-columns', {path: 'System:AppData/Chem/tests/smiles_2_columns.csv', aliases: ['smiles_2_columns'],
  description: 'two molecule columns, smiles1 and smiles2'});
dataset('activity-cliffs', {path: 'System:AppData/Chem/tests/activity_cliffs_test.csv', aliases: ['activity_cliffs_test'],
  description: '29 molecules in "smiles" with an "Activity" column'});
dataset('smiles-with-activity', {path: 'System:AppData/Chem/tests/smiles_1K_with_activities.csv', aliases: ['smiles_1K_with_activities'],
  description: '1000 molecules in "smiles" with an integer "Activity"'});
dataset('chembl-scaffolds', {path: 'System:AppData/Chem/chembl-scaffolds.csv',
  description: 'molecules with a scaffold column'});
dataset('SMILES_highlighted', {path: 'System:AppData/Chem/tests/SMILES_highlighted.csv', aliases: ['smiles-highlighted'],
  description: 'a molecule column with a scaffold to highlight'});
dataset('chem_standards', {path: 'System:AppData/Chem/chem_standards.csv', aliases: ['chem-standards'],
  description: 'salts and parent molecules for curation'});
dataset('smiles-only', {path: 'System:DemoFiles/chem/smiles_only.csv', aliases: ['smiles_only'],
  description: '1000 molecules in canonical_smiles and no other column'});
dataset('drugs-props-train', {path: 'System:AppData/Eda/drugs-props-train.csv', aliases: ['drugs_props_train'],
  description: '663 drugs with a boolean CNS column and a dozen numeric properties: the pMPO training set of the EDA package'});
dataset('ex-smarts', {path: 'System:AppData/Chem/enumerations/ex_smarts.csv', aliases: ['ex_smarts'],
  description: 'a SMARTS column of substructure patterns'});

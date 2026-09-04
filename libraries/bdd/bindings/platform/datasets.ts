/* Dataset aliases → platform locations. `{dataset}` also accepts a literal `System:…` path. */
import {dataset} from '../../src/registry.js';

dataset('spgi', {path: 'System:AppData/Chem/tests/spgi-100.csv', aliases: ['spgi-100'],
  description: 'SMILES + numeric activity, 100 rows'});
dataset('demog', {path: 'System:DemoFiles/demog.csv', description: 'the demographics demo table'});
dataset('cars', {path: 'System:DemoFiles/cars.csv'});

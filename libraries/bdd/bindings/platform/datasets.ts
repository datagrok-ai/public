/* Dataset aliases → platform locations. `{dataset}` also accepts a literal `System:…` path. */
import {dataset} from '../../src/registry.js';

dataset('spgi', {path: 'System:AppData/Chem/tests/spgi-100.csv', aliases: ['spgi-100'],
  description: 'SMILES + numeric activity, 100 rows'});
dataset('demog', {path: 'System:DemoFiles/demog.csv', description: 'the demographics demo table'});
dataset('demog-1000', {path: 'System:DemoFiles/demog-1000.csv',
  description: 'a stratified 1000-row subset of demog (same SEX / RACE / DIS_POP proportions) — the table for viewer features: every paint costs one marker per row'});
dataset('cars', {path: 'System:DemoFiles/cars.csv'});
dataset('earthquakes', {path: 'System:DemoFiles/geo/earthquakes.csv',
  description: '2426 quakes with Latitude / Longitude / Depth / Magnitude — the geo table the map viewers bind to'});
dataset('beer', {path: 'System:DemoFiles/beer.csv', description: '118 beers, 33 columns; Aroma is a long-text column, so its default filter is a text filter'});
dataset('curves', {path: 'System:DemoFiles/curves.csv', description: 'fit curves ("multiple prefit" carries the fit semantic type) next to a smiles column'});
dataset('spgi-linked1', {path: 'System:AppData/ApiTests/datasets/SPGI-linked1.csv',
  description: 'the table linked to spgi-100 by Id / Concept Id (the ApiTests package must be published)'});
dataset('spgi-linked2', {path: 'System:AppData/ApiTests/datasets/SPGI-linked2.csv',
  description: 'the table linked to spgi-linked1 by four key columns (the ApiTests package must be published)'});

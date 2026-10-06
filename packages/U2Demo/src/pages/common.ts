import {ReadonlySignal, divH, span} from '@datagrok-libraries/u2';

export const CITIES = [
  'Amsterdam', 'Athens', 'Barcelona', 'Berlin', 'Boston', 'Brussels', 'Budapest', 'Chicago',
  'Copenhagen', 'Dublin', 'Geneva', 'Helsinki', 'Kyiv', 'Lisbon', 'London', 'Madrid', 'Milan',
  'Munich', 'New York', 'Oslo', 'Paris', 'Prague', 'Rome', 'Seattle', 'Stockholm', 'Tokyo',
  'Toronto', 'Vienna', 'Warsaw', 'Zurich',
];

export const DRUGS = [
  {name: 'Aspirin', smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O'},
  {name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'},
  {name: 'Ibuprofen', smiles: 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O'},
  {name: 'Paracetamol', smiles: 'CC(=O)NC1=CC=C(C=C1)O'},
  {name: 'Naproxen', smiles: 'CC(C1=CC2=CC=C(C=C2C=C1)OC)C(=O)O'},
];

export function readout(name: string, source: string | ReadonlySignal<unknown>): HTMLElement {
  return divH([span(`${name} = `), span(source)], 'u2demo-status');
}

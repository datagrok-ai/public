import { Track } from '../../track';
import { VirtualScreeningTutorial } from './tutorials/virtual-screening';


export const tutorials = [
  VirtualScreeningTutorial,
];

export const chem = new Track(
  'Cheminformatics',
  tutorials.map((t) => new t()),
  'https://datagrok.ai/help/domains/chem/cheminformatics'
);

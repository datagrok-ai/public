import { Track } from '@datagrok-libraries/tutorials/src/track';
import { PeptidesSarTutorial } from './tutorials/peptides-sar';


export const tutorials = [
  // VirtualScreeningTutorial,
  PeptidesSarTutorial
];

export const bio = new Track(
  'Bioinformatics',
  tutorials.map((t) => new t()),
  'https://datagrok.ai/help/datagrok/solutions/domains/bio/'
);

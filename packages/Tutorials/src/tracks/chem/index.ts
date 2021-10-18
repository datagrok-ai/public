import { Track } from '../../track';
import { VirtualScreeningTutorial } from './tutorials/virtual-screening';
import { DescriptorsTutorial } from './tutorials/descriptors';


export const tutorials = [
  VirtualScreeningTutorial,
  DescriptorsTutorial,
];

export const chem = new Track(
  'Cheminformatics',
  tutorials.map((t) => new t()),
  'https://datagrok.ai/help/domains/chem/cheminformatics'
);

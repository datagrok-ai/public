import { Track } from '../../track';
import { DescriptorsTutorial } from './tutorials/descriptors';


export const tutorials = [
  DescriptorsTutorial,
];

export const chem = new Track('Cheminformatics', tutorials.map((t) => new t()), 'https://datagrok.ai/help/domains/chem/cheminformatics');

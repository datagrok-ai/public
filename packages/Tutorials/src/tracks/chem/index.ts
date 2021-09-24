import { Track } from '../../tutorial';
import { DescriptorsTutorial } from './tutorials/descriptors';


export const tutorials = [
  DescriptorsTutorial,
];

export const chem = new Track('Cheminformatics', ...tutorials.map((t) => new t()));

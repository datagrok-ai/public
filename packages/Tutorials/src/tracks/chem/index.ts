import { Track } from '../../track';
import { ActivityPredictionTutorial } from './tutorials/activity-prediction';
import { DescriptorsTutorial } from './tutorials/descriptors';


export const tutorials = [
  ActivityPredictionTutorial,
  DescriptorsTutorial,
];

export const chem = new Track(
  'Cheminformatics',
  tutorials.map((t) => new t()),
  'https://datagrok.ai/help/domains/chem/cheminformatics'
);

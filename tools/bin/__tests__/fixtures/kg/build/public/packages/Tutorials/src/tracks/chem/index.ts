import {Track} from '@datagrok-libraries/tutorials/src/track';
import {ActivityCliffsTutorial} from './tutorials/activity-cliffs';

export const tutorials = [
  ActivityCliffsTutorial,
];

export const chem = new Track(
  'Cheminformatics',
  tutorials.map((t) => new t()),
  'https://datagrok.ai/help/datagrok/solutions/domains/chem/'
);

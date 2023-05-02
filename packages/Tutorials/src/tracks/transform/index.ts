import { Track } from '@datagrok-libraries/tutorials/src/track';
import { CalculatedColumnsTutorial } from './tutorials/add-new-column';


export const tutorials = [
  CalculatedColumnsTutorial,
];

export const dataTransformation = new Track('Data transformation',
  tutorials.map((t) => new t()), '');

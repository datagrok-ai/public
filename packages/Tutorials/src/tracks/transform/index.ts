import { Track } from '@datagrok-libraries/tutorials/src/track';
import { AggregationTutorial } from './tutorials/data-aggregation';
import { CalculatedColumnsTutorial } from './tutorials/calculated-columns';


export const tutorials = [
  AggregationTutorial,
  CalculatedColumnsTutorial,
];

export const dataTransformation = new Track('Data transformation',
  tutorials.map((t) => new t()), 'https://datagrok.ai/help/transform/');

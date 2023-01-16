import { Track } from '@datagrok-libraries/tutorials/src/track';
import { DataConnectorsTutorial } from './tutorials/data-connectors';


export const tutorials = [
  DataConnectorsTutorial,
];

export const da = new Track('Data Access',
  tutorials.map((t) => new t()),
  'https://datagrok.ai/help/develop/how-to/access-data');

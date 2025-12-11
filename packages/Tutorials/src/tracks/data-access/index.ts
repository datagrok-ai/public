import { Track } from '@datagrok-libraries/tutorials/src/track';
import { DataConnectorsTutorial } from './tutorials/data-connectors';
import { StickyMetaTutorial } from './tutorials/sticky-meta';


export const tutorials = [
  DataConnectorsTutorial,
  StickyMetaTutorial,
];

export const da = new Track('Data Access',
  tutorials.map((t) => new t()),
  'https://datagrok.ai/help/develop/how-to/access-data');

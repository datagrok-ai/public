import { Track } from '@datagrok-libraries/tutorials/src/track';
import { VirtualScreeningTutorial } from './tutorials/virtual-screening';
import {RGroupsAnalysisTutorial} from './tutorials/r-groups-analysis';
import {ActivityCliffsTutorial} from './tutorials/activity-cliffs';
import {SimilarityDiversitySearchTutorial} from './tutorials/similarity-diversity-search';


export const tutorials = [
  VirtualScreeningTutorial,
  RGroupsAnalysisTutorial,
  ActivityCliffsTutorial,
  SimilarityDiversitySearchTutorial
];

export const chem = new Track(
  'Cheminformatics',
  tutorials.map((t) => new t()),
  'https://datagrok.ai/help/domains/chem/cheminformatics'
);

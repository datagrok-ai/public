import { Track } from '@datagrok-libraries/tutorials/src/track';
import { VirtualScreeningTutorial } from './tutorials/virtual-screening';
import {RGroupsAnalysisTutorial} from './tutorials/r-groups-analysis';
import {ActivityCliffsTutorial} from './tutorials/activity-cliffs';
import {SimilarityDiversitySearchTutorial} from './tutorials/similarity-diversity-search';
import {SubstructureSearchFilteringTutorial} from './tutorials/substructure-search-filtering';


export const tutorials = [
  // VirtualScreeningTutorial,
  RGroupsAnalysisTutorial,
  ActivityCliffsTutorial,
  SimilarityDiversitySearchTutorial,
  SubstructureSearchFilteringTutorial
];

export const chem = new Track(
  'Cheminformatics',
  tutorials.map((t) => new t()),
  'https://datagrok.ai/help/datagrok/solutions/domains/chem/'
);

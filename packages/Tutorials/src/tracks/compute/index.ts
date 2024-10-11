import { Track } from '@datagrok-libraries/tutorials/src/track';
import { FittingTutorial } from './tutorials/fitting-tutorial';
import { SensitivityAnalysisTutorial } from './tutorials/sensitivity-analysis-tutorial';


export const tutorials = [
  FittingTutorial,
  SensitivityAnalysisTutorial,
];

export const scientificComputing = new Track('Scientific computing', tutorials.map((t) => new t()), 'https://datagrok.ai/help/compute/');

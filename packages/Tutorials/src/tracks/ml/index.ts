import { Track } from '@datagrok-libraries/tutorials/src/track';
import { MultivariateAnalysisTutorial } from './tutorials/multivariate-analysis';
import { PredictiveModelingTutorial } from './tutorials/predictive-modeling';
import { ScriptingTutorial } from './tutorials/scripting';


export const tutorials = [
  MultivariateAnalysisTutorial,
  PredictiveModelingTutorial,
  ScriptingTutorial,
];

export const ml = new Track('Machine Learning', tutorials.map((t) => new t()), 'https://datagrok.ai/help/learn/data-science');

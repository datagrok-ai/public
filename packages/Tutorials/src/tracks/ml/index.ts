import { Track } from '../../tutorial';
import { PredictiveModelingTutorial } from './tutorials/predictive-modeling';
import { ScriptingTutorial } from './tutorials/scripting';


export const tutorials = [
  PredictiveModelingTutorial,
  ScriptingTutorial,
];

export const ml = new Track('Machine Learning', ...tutorials.map((t) => new t()));

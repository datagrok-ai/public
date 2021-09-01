import { Track } from '../../tutorial';
import { FiltersTutorial } from './tutorials/filters';
import { MultivariateAnalysisTutorial } from './tutorials/multivariate-analysis';
import { ScatterPlotTutorial } from './tutorials/scatter-plot';


export const tutorials = [
  FiltersTutorial,
  MultivariateAnalysisTutorial,
  ScatterPlotTutorial,
];
export const eda = new Track('Exploratory Data Analysis',
  ...tutorials.map((t) => new t()));

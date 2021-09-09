import { Track } from '../../tutorial';
import { FiltersTutorial } from './tutorials/filters';
import { MultivariateAnalysisTutorial } from './tutorials/multivariate-analysis';
import { ScatterPlotTutorial } from './tutorials/scatter-plot';
import { ViewersTutorial } from './tutorials/viewers-basics';


export const tutorials = [
  FiltersTutorial,
  MultivariateAnalysisTutorial,
  ScatterPlotTutorial,
  ViewersTutorial,
];
export const eda = new Track('Exploratory Data Analysis',
  ...tutorials.map((t) => new t()));

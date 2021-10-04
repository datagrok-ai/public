import { Track } from '../../track';
import { FiltersTutorial } from './tutorials/filters';
import { ScatterPlotTutorial } from './tutorials/scatter-plot';
import { ViewersTutorial } from './tutorials/viewers-basics';


export const tutorials = [
  FiltersTutorial,
  ScatterPlotTutorial,
  ViewersTutorial,
];
export const eda = new Track('Exploratory Data Analysis',
  tutorials.map((t) => new t()), 'https://datagrok.ai/help/explore/exploratory-data-analysis');

import { Track } from '../../tutorial';
import { ScatterPlotTutorial } from './tutorials/scatter-plot';


export const tutorials = [
  ScatterPlotTutorial,
];
export const eda = new Track('Exploratory Data Analysis',
  ...tutorials.map((t) => new t()));

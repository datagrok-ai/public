import { Track } from '@datagrok-libraries/tutorials/src/track';
import { DashboardTutorial } from './tutorials/dashboard';
import { FiltersTutorial } from './tutorials/filters';
import { ScatterPlotTutorial } from './tutorials/scatter-plot';
import { ViewersTutorial } from './tutorials/viewers-basics';
import { EmbeddedViewersTutorial } from './tutorials/embedded-viewers';
import { GridTutorial } from './tutorials/grid';


export const tutorials = [
  ScatterPlotTutorial,
  GridTutorial,
  ViewersTutorial,
  FiltersTutorial,
  DashboardTutorial,
  EmbeddedViewersTutorial,
];

export const eda = new Track('Exploratory Data Analysis',
  tutorials.map((t) => new t()), 'https://datagrok.ai/help/explore/exploratory-data-analysis');

import * as DG from 'datagrok-api/dg';
import {Property} from "datagrok-api/src/entities";
import {TYPE} from "datagrok-api/src/const";

export const FIT_TYPE = 'fit';

export enum ChartPointStyle {
  points = 'points',
  boxplot = 'boxplot',
  hide = 'hide'
}

/** Chart options. For fitted curves, this object is stored in the grid column tags and is
 * used by the renderer. */
export interface IChartOptions {
  minX? : number;
  minY? : number;
  maxX? : number;
  maxY? : number;

  xAxisName?: string;
  yAxisName?: string;

  xAxisColumnName?: string;
  yAxisColumnName?: string;

  /** If true, clicking on the point toggles its outlier status and causes curve refitting */
  clickToToggle?: boolean;

  autoFit?: boolean;
  fitFunction?: string;
  showFitLine?: boolean;
  showPoints?: boolean;
  pointStyle?: ChartPointStyle;
  color?: string;
}

export const fitColumnProperties: Property[] = [
  // Style and zoom
  Property.js('minX', TYPE.FLOAT, {description: 'Minimum value of the X axis', nullable: true}),
  Property.js('minY', TYPE.FLOAT, {description: 'Minimum value of the Y axis', nullable: true}),
  Property.js('maxX', TYPE.FLOAT, {description: 'Maximum value of the X axis', nullable: true}),
  Property.js('maxY', TYPE.FLOAT, {description: 'Maximum value of the Y axis', nullable: true}),
  Property.js('xAxisName', TYPE.STRING, {description: 'Label to show on the X axis. If not specified, corresponding data column name is used', nullable: true}),
  Property.js('yAxisName', TYPE.STRING, {description: 'Label to show on the Y axis. If not specified, corresponding data column name is used', nullable: true}),

  // Fitting
  Property.js('clickToToggle', TYPE.BOOL, {category: 'Fitting', description: 'If true, clicking on the point toggles its outlier status and causes curve refitting', nullable: true}),
  Property.js('autoFit', TYPE.BOOL, {category: 'Fitting', description: 'Perform fitting on-the-fly', defaultValue: true}),
  Property.js('showFitLine', TYPE.BOOL, {category: 'Fitting', description: 'Whether the fit line should be rendered', defaultValue: true}),
  Property.js('showPoints', TYPE.BOOL, {category: 'Fitting', description: 'Whether points should be rendered', defaultValue: true}),
  Property.js('fitFunction', TYPE.STRING, {category: 'Fitting', choices: ['sigmoid', 'linear'], defaultValue: 'sigmoid'}),

  // Rendering
  Property.js('pointStyle', TYPE.STRING, {category: 'Rendering', choices: [ChartPointStyle.points, ChartPointStyle.boxplot], defaultValue: ChartPointStyle.points}),
  Property.js('color', TYPE.STRING, {category: 'Rendering', defaultValue: DG.Color.toHtml(DG.Color.scatterPlotMarker), nullable: true}),
];


export const fitResultProperties: Property[] = [
  Property.js('sigma', TYPE.FLOAT),
  Property.js('rSquared', TYPE.FLOAT),
  Property.js('hillSlope', TYPE.FLOAT),
  Property.js('classification', TYPE.STRING),
];


/** Creates new object with the default values specified in the properties */
function createFromProperties(props: Property[]): any {
  const o: any = {};
  for (let p of props)
    if (p.defaultValue != null)
      o[p.name] = p.defaultValue;
  return o;
}


export function getColumnFitOptions(gridColumn: DG.GridColumn): IChartOptions {
  return gridColumn.tags['fitOptions'] ??= createFromProperties(fitColumnProperties);
}
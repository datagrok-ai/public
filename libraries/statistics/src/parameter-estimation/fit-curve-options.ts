import {Property} from "datagrok-api/src/entities";
import {TYPE} from "datagrok-api/src/const";

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

  pointStyle?: ChartPointStyle;
}

export const optionsProperties: Property[] = [
  Property.js('minX', TYPE.FLOAT, {description: 'Minimum value of the X axis', nullable: true}),
  Property.js('minY', TYPE.FLOAT, {description: 'Minimum value of the Y axis', nullable: true}),
  Property.js('maxX', TYPE.FLOAT, {description: 'Maximum value of the X axis', nullable: true}),
  Property.js('maxY', TYPE.FLOAT, {description: 'Maximum value of the Y axis', nullable: true}),
  Property.js('xAxisName', TYPE.STRING, {description: 'Label to show on the X axis. If not specified, corresponding data column name is used', nullable: true}),
  Property.js('yAxisName', TYPE.STRING, {description: 'Label to show on the Y axis. If not specified, corresponding data column name is used', nullable: true}),
  Property.js('clickToToggle', TYPE.BOOL, {category: 'Behavior', description: 'If true, clicking on the point toggles its outlier status and causes curve refitting', nullable: true}),
];

export interface IChartSeries {

}

export interface IChartData {

}
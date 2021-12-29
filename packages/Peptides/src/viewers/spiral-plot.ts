import * as DG from 'datagrok-api/dg';
import {_toJson} from 'datagrok-api/src/utils';

import {assert, argSort} from '@datagrok-libraries/utils/src/operations';
import {Options} from '@datagrok-libraries/utils/src/type-declarations';

const api = <any>window;

/**
 * Draws 2D scatter plot from 1D series.
 *
 * @export
 * @class SpiralPlot
 * @extends {DG.ScatterPlotViewer}
 */
export class SpiralPlot extends DG.ScatterPlotViewer {
    static axesNames = ['~X', '~Y'];
    static valuesKey = 'valuesColumnName';

    /**
     * Calculates coordinates of the projection into a spiral.
     *
     * @static
     * @param {DG.DataFrame} t Source data frame.
     * @param {Options} options Options to read values column name from. Must include {valuesColumnName: string}.
     * @return {DG.DataFrame} Updated dataframe.
     * @memberof SpiralPlot
     */
    static updateCoordinates(t: DG.DataFrame, options: Options): DG.DataFrame {
      assert(options[SpiralPlot.valuesKey] != undefined);

      const values = t.getCol(options[SpiralPlot.valuesKey]).getRawData() as Float32Array;
      const columns = _calcSpiralProjection(values);
      const cdf = DG.DataFrame.fromColumns(
        Array.from(columns).map((v, i) => DG.Column.fromFloat32Array(SpiralPlot.axesNames[i], v)),
      );
      return _updateCoordinates(t, cdf);
    }

    /**
     * Creates new SpiralPlot from a data frame with selected values column.
     *
     * @static
     * @param {DG.DataFrame} t A data frame.
     * @param {Options} options Controlling options.
     * @return {SpiralPlot} The plot.
     * @memberof SpiralPlot
     */
    static fromTable(t: DG.DataFrame, options: Options): SpiralPlot {
      t = SpiralPlot.updateCoordinates(t, options);
      [options.x, options.y] = SpiralPlot.axesNames;
      options.color = options[SpiralPlot.valuesKey];
      return new SpiralPlot(api.grok_Viewer_ScatterPlot(t.dart, _toJson(options)));
    }
}

/**
 * Calculates 2D projection of 1D series as a spiral.
 *
 * @param {(number[] | Float32Array)} values The series.
 * @return {[Float32Array, Float32Array]} X and Y componenets of the projection.
 */
function _calcSpiralProjection(values: number[] | Float32Array): [Float32Array, Float32Array] {
  const nItems = values.length;
  const order = argSort(Array.from(values), true);
  const maxV = values[order[0]];
  const X = new Float32Array(nItems).fill(0);
  const Y = new Float32Array(nItems).fill(0);

  for (const i of order) {
    const v = maxV - values[i];
    X[i] = v * Math.cos(Math.PI * v) - Math.random() * 1.5 + 0.75;
    Y[i] = v * Math.sin(Math.PI * v) - Math.random() * 1.5 + 0.75;
  }
  return [X, Y];
}

/**
 * Adds new columns from one data frame into another one.
 *
 * @param {DG.DataFrame} table Destination data frame.
 * @param {DG.DataFrame} coords Source data frame.
 * @return {DG.DataFrame} Updated data frame.
 */
function _updateCoordinates(table: DG.DataFrame, coords: DG.DataFrame): DG.DataFrame {
  const coordsColNames: string[] = coords.columns.names();
  const tableColNames: string[] = table.columns.names();
  const restColNames = tableColNames.filter((v: string) => !coordsColNames.includes(v));

  if (tableColNames.length == restColNames.length) {
    for (const col of coords.columns) {
      table.columns.add(col);
    }
    return table;
  }
  return table.join(coords, coordsColNames, coordsColNames, restColNames, [], 'right', false);
}

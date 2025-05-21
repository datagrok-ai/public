/* eslint-disable max-len */
// Visualizer of results of multi-objective optimization

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {isDimValid} from './utils';
import {OPT_TYPE} from './defs';
import {getOutputPalette, INPUT_COLOR_PALETTE} from './ui-tools';

function singleInputSingleOutputViewer(table: DG.DataFrame): DG.Viewer[] {
  const colsCount = table.columns.length;

  if (colsCount !== 2)
    throw new Error(`Failed to create "single-input-single-output" viewer. Incorrect columns count: ${colsCount}, expected: 2`);

  return [DG.Viewer.lineChart(table)];
}

function singleInputDoubleOutputViewer(table: DG.DataFrame): DG.Viewer[] {
  const colsCount = table.columns.length;

  if (colsCount <= 3)
    throw new Error(`Failed to create "single-input-dounle-output" viewer. Incorrect columns count: ${colsCount}, expected: 2`);

  return [DG.Viewer.lineChart(table)];
}

export function visualize(view: DG.TableView, inputDim: number, outputDim: number): void {
  if (!isDimValid(inputDim))
    throw new Error(`Invalid inputs dimension: ${inputDim}. The optimizer issue`);

  if (!isDimValid(outputDim))
    throw new Error(`Invalid outputs dimension: ${outputDim}. The optimizer issue`);

  const table = view.dataFrame;

  //if ((inputDim === 1) && (outputDim === 1))
  //return singleInputSingleOutputViewer(table);
}

export class Visualizer {
  private view: DG.TableView;
  private table: DG.DataFrame;
  private inputDim: number;
  private outputDim: number;
  private type: OPT_TYPE;

  constructor(view: DG.TableView, inputDim: number, outputDim: number, type: OPT_TYPE) {
    if (!isDimValid(inputDim))
      throw new Error(`Invalid inputs dimension: ${inputDim}. The optimizer issue`);

    if (!isDimValid(outputDim))
      throw new Error(`Invalid outputs dimension: ${outputDim}. The optimizer issue`);

    const table = view.table;

    if (table === null)
      throw new Error('Null dataframe with optimization resuls. The optimizer issue');

    this.view = view;
    this.table = table;
    this.inputDim = inputDim;
    this.outputDim = outputDim;
    this.type = type;
  }

  public visualize(): void {
    this.setColorCoding();
  }

  private setColorCoding(): void {
    const cols = this.table.columns;

    for (let i = 0; i < this.inputDim; ++i) {
      const col = cols.byIndex(i);
      col.meta.colors.setLinear(
        INPUT_COLOR_PALETTE,
        {min: col.stats.min, max: col.stats.max},
      );
    }

    for (let i = this.inputDim; i < this.inputDim + this.outputDim; ++i) {
      const col = cols.byIndex(i);
      col.meta.colors.setLinear(
        getOutputPalette(this.type),
        {min: col.stats.min, max: col.stats.max},
      );
    }
  }
}

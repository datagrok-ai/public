/* eslint-disable max-len */
// Visualizer of results of multi-objective optimization

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {isDimValid} from './utils';
import {OPT_TYPE} from './defs';
import {DOCK_RATIO, getOutputPalette, INPUT_COLOR_PALETTE, TITLE} from './ui-tools';

export class Visualizer {
  private view: DG.TableView;
  private table: DG.DataFrame;
  private colNames: string[];
  private inputDim: number;
  private outputDim: number;
  private type: OPT_TYPE;
  private gridNode: DG.DockNode;

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
    this.colNames = table.columns.names();

    const node = this.view.dockManager.findNode(this.view.grid.root);
    if (node === undefined)
      throw new Error('The view is corrupted: not found the grid');

    this.gridNode = node;
  }

  public visualize(): DG.Viewer[] {
    this.setColorCoding();
    this.sortResults();

    return this.addViewers();
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

  private sortResults(): void {
    if (this.outputDim !== 1)
      return;

    const outColName = this.table.columns.byIndex(this.inputDim).name;
    this.view.grid.sort([outColName], [this.type === OPT_TYPE.MIN]);
  }

  private addViewers(): DG.Viewer[] {
    // Multiple objectives case
    if (this.outputDim > 2)
      return this.addMultipleInputsMultipleOutputViewers();

    // Double objectives case
    if (this.outputDim === 2)
      return this.addMultipleInputsMultipleOutputViewers();

    // Single objective case
    if (this.inputDim === 1)
      return this.addSingleInputSingleOutputViewer();

    return this.addMultipleInputsSingeOutputViewers();
  }

  private getScatter(): DG.Viewer {
    return DG.Viewer.scatterPlot(this.table, {
      xColumnName: this.colNames[0],
      yColumnName: this.colNames[1],
      colorColumnName: this.colNames[this.inputDim],
    });
  }

  private getParetoScatter(): DG.Viewer {
    return DG.Viewer.scatterPlot(this.table, {
      title: TITLE.PARETO_FRONT,
      xColumnName: this.colNames[this.inputDim],
      yColumnName: this.colNames[this.inputDim + 1],
      colorColumnName: (this.inputDim === 1) ? this.colNames[0] : (this.outputDim >= 3 ? this.colNames[this.inputDim + 2] : undefined),
      labelColumnNames: (this.inputDim === 2) ? [this.colNames[0], this.colNames[1]] : undefined,
    });
  }

  private getPc(): DG.Viewer {
    return DG.Viewer.pcPlot(this.table, {showColorSelector: false});
  }

  private addSingleInputSingleOutputViewer(): DG.Viewer[] {
    const chart = this.getScatter();

    this.view.dockManager.dock(
      chart,
      DG.DOCK_TYPE.LEFT,
      this.gridNode,
      undefined,
      DOCK_RATIO.SINGLE_SINGLE_SCATTER,
    );

    return [chart];
  }

  private addMultipleInputsSingeOutputViewers(): DG.Viewer[] {
    const scatter = this.getScatter();
    const scatterNode = this.view.dockManager.dock(
      scatter,
      DG.DOCK_TYPE.LEFT,
      this.gridNode,
      undefined,
      DOCK_RATIO.SINGLE_SINGLE_SCATTER,
    );

    const pcPlot = this.getPc();
    this.view.dockManager.dock(
      pcPlot,
      DG.DOCK_TYPE.DOWN,
      scatterNode,
      undefined,
      DOCK_RATIO.MULT_SINGLE_PC,
    );

    return [scatter, pcPlot];
  }

  private addMultipleInputsMultipleOutputViewers(): DG.Viewer[] {
    const scatter = this.getParetoScatter();
    const scatterNode = this.view.dockManager.dock(
      scatter,
      DG.DOCK_TYPE.LEFT,
      this.gridNode,
      undefined,
      DOCK_RATIO.SINGLE_SINGLE_SCATTER,
    );

    const pcPlot = this.getPc();
    this.view.dockManager.dock(
      pcPlot,
      DG.DOCK_TYPE.DOWN,
      scatterNode,
      undefined,
      DOCK_RATIO.MULT_SINGLE_PC,
    );

    return [scatter, pcPlot];
  }
}

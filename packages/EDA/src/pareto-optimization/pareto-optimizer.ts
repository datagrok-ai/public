import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/pareto.css';
import {NumericFeature, OPT_TYPE, NumericArray, DIFFERENCE, RATIO, COL_NAME, PC_MAX_COLS, ColorOpt} from './defs';
import {getColorScaleDiv, getOutputPalette} from './utils';

/** Pareto front optimization app */
export class ParetoOptimizer {
  private df: DG.DataFrame;
  private numCols: DG.Column[];
  private numColNames: string[] = [];
  private numColsCount: number;
  private rowCount: number;
  private features = new Map<string, NumericFeature>();
  private view: DG.TableView;
  private pcPlot: DG.Viewer<DG.IPcPlotSettings>;
  private toUpdatePcCols = false;

  private paretoFrontViewer: DG.Viewer;
  private resultColName: string;

  private intervalId: NodeJS.Timeout | null = null;

  private inputsMap = new Map<string, DG.InputBase>();

  private subs: any[];
  private pcPlotNode: DG.DockNode | null = null;
  private inputFormNode: DG.DockNode | null = null;

  constructor(df: DG.DataFrame) {
    this.df = df;
    const cols = df.columns;
    const colList = cols.toList();
    this.numCols = colList.filter((col) => col.isNumerical);
    this.numColNames = this.numCols.map((col) => col.name);
    this.numColsCount = this.numCols.length;
    this.rowCount = df.rowCount;
    this.view = grok.shell.getTableView(df.name);

    this.paretoFrontViewer = DG.Viewer.fromType('Pareto front', df);

    const paretoFrontViewerNode = this.view.dockManager.dock(
      this.paretoFrontViewer,
      DG.DOCK_TYPE.RIGHT,
      null,
      undefined,
      RATIO.VIEWER,
    );

    this.pcPlot = DG.Viewer.pcPlot(df, {legendPosition: 'Top'});
    const gridNode = this.view.dockManager.findNode(this.view.grid.root) ?? paretoFrontViewerNode;
    this.pcPlotNode = this.view.dockManager.dock(this.pcPlot, DG.DOCK_TYPE.DOWN, gridNode, undefined, RATIO.VIEWER);
    this.toUpdatePcCols = this.numColNames.length > PC_MAX_COLS;

    this.resultColName = this.df.columns.getUnusedName(COL_NAME.OPT);
    this.showResultOptCol();
    this.subs = this.getSubscriptions();
  } // constructor

  private isApplicable(): boolean {
    if (this.rowCount < 1) {
      grok.shell.warning('Cannot compute Pareto front: the table is empty.');
      return false;
    }

    if (this.numColsCount < 2) {
      grok.shell.warning('Cannot compute Pareto front: at least two numeric columns are required.');
      return false;
    }

    return true;
  } // isApplicable

  public run(): void {
    if (!this.isApplicable())
      return;

    this.buildInputsForm();
    this.computeParetoFront();
    this.updateVisualization();
  } // run

  private getSubscriptions() {
    const subs = [
      this.paretoFrontViewer.onDetached.subscribe(() => {
        if (this.pcPlotNode !== null) {
          this.view.dockManager.close(this.pcPlotNode);
          this.pcPlotNode = null;
        }

        if (this.inputFormNode !== null) {
          this.view.dockManager.close(this.inputFormNode);
          this.inputFormNode = null;
        }

        this.numCols.forEach((col) => col.colors.setDisabled());
        this.features.clear();
      }),
    ];

    return subs;
  } // getSubscriptions

  private buildInputsForm(): void {
    const form = ui.form([]);
    form.classList.add('pareto-input-form');
    form.append(ui.h1('Optimize'));

    this.numCols.forEach((col, idx) => {
      const feature: NumericFeature = {
        toOptimize: this.numColsCount - idx - 1 < DIFFERENCE,
        optType: OPT_TYPE.MIN,
      };

      const name = col.name;

      const optimizationTypeInput = ui.input.choice(name, {
        value: feature.toOptimize ? feature.optType : null,
        nullable: true,
        items: [null, OPT_TYPE.MIN, OPT_TYPE.MAX],
        onValueChanged: (val) => {
          if (val == null)
            feature.toOptimize = false;
          else {
            feature.toOptimize = true;
            feature.optType = val;
          }

          this.computeParetoFront();
          this.updateVisualization();
        },
      });
      ui.tooltip.bind(optimizationTypeInput.input, () => {
        if (feature.toOptimize)
          return ui.markdown(`M${feature.optType.slice(1)} **${name}** during Pareto optimization`);

        return ui.markdown(`Ignore **${name}** during Pareto optimization`);
      });

      this.inputsMap.set(name, optimizationTypeInput);

      form.append(optimizationTypeInput.root);
      this.features.set(name, feature);
    });

    this.inputFormNode = this.view.dockManager.dock(form, DG.DOCK_TYPE.LEFT, null, undefined, RATIO.FORM);
  } // buildInputsForm

  private computeParetoFront(): void {
    const minimizeColumnNames: string[] = [];
    const maximizeColumnNames: string[] = [];

    this.features.forEach((fea, name) => {
      if (fea.toOptimize) {
        if (fea.optType === OPT_TYPE.MIN)
          minimizeColumnNames.push(name);
        else
          maximizeColumnNames.push(name);
      }
    });

    this.paretoFrontViewer.setOptions({
      minimizeColumnNames: minimizeColumnNames,
      maximizeColumnNames: maximizeColumnNames,
    });
  } // computeParetoFront

  private updatePcPlot(colNames: string[], colorOpt: ColorOpt): void {
    this.pcPlot.setOptions(colorOpt);

    // update value columns: check that optimized cols are included
    if (this.toUpdatePcCols) {
      const prevColNames = this.pcPlot.getOptions().look['columnNames'];

      let toUpdatePcPlotColNames = false;

      colNames.forEach((name) => {
        if (!prevColNames.includes(name))
          toUpdatePcPlotColNames = true;
      });

      if (toUpdatePcPlotColNames) {
        const valColNames = [...colNames];
        const notIncluded = this.numColNames.filter((name) => !valColNames.includes(name));
        valColNames.push(...notIncluded.slice(0, PC_MAX_COLS - colNames.length));
        this.pcPlot.setOptions({columnNames: valColNames});
      }
    }
  } // updatePcPlot

  private updateVisualization(): void {
    const colNames: string[] = [];

    this.features.forEach((fea, name) => {
      if (fea.toOptimize)
        colNames.push(name);
    });

    const colorOpt: ColorOpt = {'colorColumnName': (colNames.length > 0) ? this.resultColName: undefined};
    this.updatePcPlot(colNames, colorOpt);
    this.markOptColsWithColor();
    this.updateTooltips();
  } // updateVisualization

  private showResultOptCol(): void {
    // show a column with the results, once it is added
    this.intervalId = setInterval(() => {
      const gridCol = this.view.grid.columns.byName(this.resultColName);

      if (gridCol !== null) {
        gridCol.visible = true;
        this.stopChecking();
      }
    }, 1000);
  } // showResultOptCol

  private stopChecking(): void {
    if (this.intervalId) {
      clearInterval(this.intervalId);
      this.intervalId = null;
    }
  }

  private markOptColsWithColor(): void {
    this.numCols.forEach((col) => col.colors.setDisabled());

    this.features.forEach((fea, name) => {
      if (!fea.toOptimize)
        return;

      const col = this.df.col(name);

      if (col != null)
        col.colors.setLinear(getOutputPalette(fea.optType), {min: col.stats.min, max: col.stats.max});
    });
  } // markOptColsWithColor

  private updateTooltips(): void {
    const features = this.features;

    this.view.grid.onCellTooltip(function(cell, x, y) {
      if (cell.isColHeader) {
        const cellCol = cell.tableColumn;
        if (cellCol) {
          const name = cell.tableColumn.name;
          const feature = features.get(name);

          if (feature !== undefined) {
            const elems = [ui.markdown(`**${name}**`)];

            if (feature.toOptimize) {
              elems.push(ui.markdown(`This feature is **${feature.optType}d** during Pareto optimization.`));
              elems.push(getColorScaleDiv(feature.optType));
            }

            ui.tooltip.show(ui.divV(elems), x, y);

            return true;
          }

          return false;
        }
      }
    });
  } // updateTooltips
} // ParetoOptimizer

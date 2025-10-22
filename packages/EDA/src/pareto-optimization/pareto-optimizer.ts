import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/pareto.css';
import {paretoMaskFromCoordinates} from './pareto-computations';
import {NumericFeature, OPT_TYPE, NumericArray, DIFFERENCE, RATIO, COL_NAME,
  PC_MAX_COLS, AXIS_NAMES, ColorOpt, AXIS_NAMES_3D, LABEL, SIZE, SCATTER_ROW_LIM, SCATTER3D_ROW_LIM} from './defs';
import {getColorScaleDiv, getOutputPalette, PALETTE} from './utils';

export class ParetoOptimizer {
  private df: DG.DataFrame;
  private numCols: DG.Column[];
  private numColNames: string[] = [];
  private numColsCount: number;
  private rowCount: number;
  private features = new Map<string, NumericFeature>();
  private view: DG.TableView;
  private scatter: DG.ScatterPlotViewer;
  private resultColName: string;
  private sizeColName: string;
  private labelColumnNames: string[];
  private scatter3d: DG.Viewer<DG.IScatterPlot3dSettings> | null = null;
  private pcPlot: DG.Viewer<DG.IPcPlotSettings>;
  private toUpdatePcCols = false;
  private toChangeScatterMarkerSize = false;
  private toChange3dScatterMarkerSize = false;
  private toAddSizeCol = false;
  private useAxesInput = ui.input.bool('Pareto Axes', {
    value: true,
    tooltipText: 'Use optimized variables as scatter plot axes',
    onValueChanged: (val) => {
      if (val)
        this.updateVisualization();
    },
  });

  constructor(df: DG.DataFrame) {
    this.df = df;
    const cols = df.columns;
    const colList = cols.toList();
    this.numCols = colList.filter((col) => col.isNumerical);
    this.numColNames = this.numCols.map((col) => col.name);
    this.numColsCount = this.numCols.length;
    this.rowCount = df.rowCount;
    this.view = grok.shell.getTableView(df.name);
    this.labelColumnNames = this.getLabelColNames();
    this.scatter = DG.Viewer.scatterPlot(df, {
      labelColumnNames: this.labelColumnNames,
      legendPosition: 'Top',
      markerType: DG.MARKER_TYPE.CIRCLE,
      autoLayout: false,
      showSizeSelector: false,
    });

    this.toChangeScatterMarkerSize = this.rowCount > SCATTER_ROW_LIM;
    this.toChange3dScatterMarkerSize = this.rowCount > SCATTER3D_ROW_LIM;
    this.toAddSizeCol = this.toChangeScatterMarkerSize || this.toChange3dScatterMarkerSize;

    if (this.toChangeScatterMarkerSize)
      this.scatter.setOptions({markerMinSize: SIZE.NON_OPT, markerMaxSize: SIZE.OPTIMAL});

    const scatterNode = this.view.dockManager.dock(this.scatter, DG.DOCK_TYPE.RIGHT, null, undefined, RATIO.VIEWER);

    this.pcPlot = DG.Viewer.pcPlot(df, {legendPosition: 'Top'});
    const gridNode = this.view.dockManager.findNode(this.view.grid.root) ?? scatterNode;
    this.view.dockManager.dock(this.pcPlot, DG.DOCK_TYPE.DOWN, gridNode, undefined, RATIO.VIEWER);
    this.toUpdatePcCols = this.numColNames.length > PC_MAX_COLS;

    // if (this.numColsCount > 2) {
    // // Disabled the use of 3d scatter due to the bug: https://reddata.atlassian.net/browse/GROK-19071
    //   this.scatter3d = DG.Viewer.scatterPlot3d(df);
    //   this.view.dockManager.dock(this.scatter3d, DG.DOCK_TYPE.FILL, scatterNode);
    // }

    this.resultColName = cols.getUnusedName(COL_NAME.OPT);
    this.sizeColName = cols.getUnusedName(COL_NAME.SIZE);

    const ribPanels = this.view.getRibbonPanels();
    ribPanels.push([this.useAxesInput.root]);
    this.view.setRibbonPanels(ribPanels);

    this.updateUseAxesInput();
  } // constructor

  private getLabelColNames(): string[] {
    return this.df.columns.toList().filter((col) => {
      if (col.type !== DG.COLUMN_TYPE.STRING)
        return false;

      return col.categories.length === this.rowCount;
    }).map((col) => col.name);
  } // getLabelColNames

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

      const typeInp = ui.input.choice(name, {
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
      ui.tooltip.bind(typeInp.input, () => {
        if (feature.toOptimize)
          return ui.markdown(`M${feature.optType.slice(1)} **${name}** during Pareto optimization`);

        return ui.markdown(`Ignore **${name}** during Pareto optimization`);
      });

      //typeInp.input.hidden = !feature.toOptimize;

      const enableInp = ui.input.toggle('', {
        value: feature.toOptimize,
        onValueChanged: (val) => {
          feature.toOptimize = val;
          typeInp.input.hidden = !val;
          this.computeParetoFront();
          this.updateVisualization();
        },
      });
      ui.tooltip.bind(enableInp.root, () => `${feature.toOptimize ? 'Dis' : 'En'}able optimization for "${name}"`);

      enableInp.root.classList.add('pareto-switch-input');

      //typeInp.root.insertBefore(enableInp.root, typeInp.captionLabel);

      form.append(typeInp.root);
      this.features.set(name, feature);
    });

    this.view.dockManager.dock(form, DG.DOCK_TYPE.LEFT, null, undefined, RATIO.FORM);
  } // buildInputsForm

  private computeParetoFront(): void {
    const data: NumericArray[] = [];
    const sense: OPT_TYPE[] = [];

    this.features.forEach((fea, name) => {
      if (fea.toOptimize) {
        data.push(this.df.col(name)!.getRawData());
        sense.push(fea.optType);
      }
    });

    if (data.length > 0) {
      const mask = paretoMaskFromCoordinates(data, sense, this.rowCount);
      const colOpt = DG.Column.fromStrings(this.resultColName, mask.map((res) => res ? LABEL.OPTIMAL : LABEL.NON_OPT));
      this.df.columns.remove(this.resultColName, true);
      this.df.columns.add(colOpt);
      this.markResColWithColor(colOpt);

      if (this.toAddSizeCol) {
        const sizeCol = DG.Column.fromInt32Array(
          this.sizeColName,
          new Int32Array(mask.map((res) => res ? SIZE.OPTIMAL : SIZE.NON_OPT)),
        );
        this.df.columns.remove(this.sizeColName, true);
        this.df.columns.add(sizeCol);
      }
    } else {
      this.df.columns.remove(this.resultColName, true);
      this.df.columns.remove(this.sizeColName, true);
    }
  } // computeParetoFront

  private updateScatter(plot: DG.Viewer | null, axisNames: string[], colNames: string[], colorOpt: ColorOpt): void {
    if (plot == null)
      return;

    plot.setOptions(colorOpt);

    // update axis
    if (!this.useAxesInput.value)
      return;

    const prevAxisCols = axisNames.map((axis) => plot.getOptions().look[axis]);

    let toUpdate = false;

    colNames.forEach((name) => {
      if (!prevAxisCols.includes(name))
        toUpdate = true;
    });

    if (toUpdate) {
      axisNames.forEach((axis, idx) => {
        const opt: Record<string, string> = {};
        opt[axis] = colNames[idx];
        plot.setOptions(opt);
      });
    }
  } // updateScatter

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

  private setMarkers(): void {
    if (this.toChangeScatterMarkerSize) {
      this.scatter.setOptions({
        'sizeColumnName': this.sizeColName,
        'markerMinSize': SIZE.NON_OPT,
        'markerMaxSize': SIZE.OPTIMAL,
        'markerType': DG.MARKER_TYPE.CIRCLE,
      });
    }

    if (this.toChange3dScatterMarkerSize && (this.scatter3d !== null))
      this.scatter3d.setOptions({'size': this.sizeColName});
  } // setMarkers

  private updateVisualization(): void {
    const colNames: string[] = [];

    this.features.forEach((fea, name) => {
      if (fea.toOptimize)
        colNames.push(name);
    });

    const colorOpt: ColorOpt = {'colorColumnName': (colNames.length > 0) ? this.resultColName: undefined};
    this.updateScatter(this.scatter, AXIS_NAMES, colNames, colorOpt);
    this.updateScatter(this.scatter3d, AXIS_NAMES_3D, colNames, colorOpt);
    this.updatePcPlot(colNames, colorOpt);
    this.setMarkers();
    this.hideSizeCol();
    this.markOptColsWithColor();
    this.updateTooltips();
  } // updateVisualization

  private hideSizeCol(): void {
    const gridCol = this.view.grid.columns.byName(this.sizeColName);

    if (gridCol !== null)
      gridCol.visible = false;
  } // hideSizeCol

  private markResColWithColor(col: DG.Column): void {
    col.colors.setCategorical({
      'optimal': '#2ca02c',
      'non-optimal': '#e3e3e3',
    });
  } // markResColWithColor

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

  private updateUseAxesInput(): void {
    this.useAxesInput.captionLabel.style.fontSize = '10pt';
    this.useAxesInput.captionLabel.style.position = 'relative';
    this.useAxesInput.captionLabel.style.top = '-2px';
    this.useAxesInput.root.style.paddingTop = '5px';
  } // updateUseAxesInput
} // ParetoOptimizer

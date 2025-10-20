import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/pareto.css';
import {paretoMaskFromCoordinates} from './pareto-computations';
import {NumericFeature, OPT_TYPE, NumericArray, DIFFERENCE, RATIO, OPTIMALITY_COL_NAME,
  PC_MAX_COLS, scatterAxisNames,
  ParetoLabel,
  ColorOpt} from './defs';


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
  private labelColumnNames: string[];
  private scatter3d: DG.Viewer<DG.IScatterPlot3dSettings> | null = null;
  private pcPlot: DG.Viewer<DG.IPcPlotSettings>;
  private toUpdatePcCols = false;

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
    this.scatter = DG.Viewer.scatterPlot(df, {labelColumnNames: this.labelColumnNames, legendPosition: 'Top'});
    const scatterNode = this.view.dockManager.dock(this.scatter, DG.DOCK_TYPE.RIGHT, null, undefined, RATIO.VIEWER);

    this.pcPlot = DG.Viewer.pcPlot(df, {legendPosition: 'Top'});
    const gridNode = this.view.dockManager.findNode(this.view.grid.root) ?? scatterNode;
    this.view.dockManager.dock(this.pcPlot, DG.DOCK_TYPE.DOWN, gridNode, undefined, RATIO.VIEWER);
    this.toUpdatePcCols = this.numColNames.length > PC_MAX_COLS;

    if (this.numColsCount > 2) {
      this.scatter3d = DG.Viewer.scatterPlot3d(df);
      this.view.dockManager.dock(this.scatter3d, DG.DOCK_TYPE.FILL, scatterNode);
    }

    this.resultColName = cols.getUnusedName(OPTIMALITY_COL_NAME);
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
        value: feature.optType,
        nullable: false,
        items: [OPT_TYPE.MIN, OPT_TYPE.MAX],
        onValueChanged: (val) => {
          feature.optType = val;
          this.computeParetoFront();
          this.updateVisualization();
        },
      });
      ui.tooltip.bind(typeInp.input, 'Type of optimization');

      typeInp.input.hidden = !feature.toOptimize;

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

      typeInp.root.insertBefore(enableInp.root, typeInp.captionLabel);

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
      const colOpt = DG.Column.fromStrings(this.resultColName, mask);
      this.df.columns.remove(this.resultColName, true);
      this.df.columns.add(colOpt);
    } else
      this.df.columns.remove(this.resultColName, true);
  } // computeParetoFront

  private updateScatter(colNames: string[], colorOpt: ColorOpt): void {
    this.scatter.setOptions(colorOpt);

    // update axis
    const prevAxisCols = scatterAxisNames.map((axis) => this.scatter.getOptions().look[axis]);

    let toUpdateScatterAxisCols = false;

    colNames.forEach((name) => {
      if (!prevAxisCols.includes(name))
        toUpdateScatterAxisCols = true;
    });

    if (toUpdateScatterAxisCols) {
      scatterAxisNames.forEach((axis, idx) => {
        const opt: Record<string, string> = {};
        opt[axis] = colNames[idx];
        this.scatter.setOptions(opt);
      });
    }
  } // updateScatter

  private update3dScatter(colNames: string[], colorOpt: ColorOpt): void {
    if (this.scatter3d !== null) {
      this.scatter3d.setOptions(colorOpt);

      if (colNames.length > 0)
        this.scatter3d.setOptions({xColumnName: colNames[0]});

      if (colNames.length > 1)
        this.scatter3d.setOptions({yColumnName: colNames[1]});

      if (colNames.length > 2)
        this.scatter3d.setOptions({zColumnName: colNames[2]});
    }
  } // update3dScatter

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

    this.updateScatter(colNames, colorOpt);
    this.update3dScatter(colNames, colorOpt);
    this.updatePcPlot(colNames, colorOpt);
  } // updateVisualization
} // ParetoOptimizer

export function runParetoOptimizer(): void {
  const df: DG.DataFrame | null = grok.shell.t;
  if (df === null) {
    grok.shell.warning('No dataframe is opened');
    return;
  }

  const cols = df.columns.toList();
  const numCols = cols.filter((col) => col.isNumerical);
  if (numCols.length < 1) {
    grok.shell.warning('No numeric columns to be optimized');
    return;
  }

  const numColsCount = numCols.length;
  const rowCount = df.rowCount;

  const features = new Map<string, NumericFeature>();
  const form = ui.form([]);
  form.classList.add('pareto-input-form');

  //let optCol = df.col('__Optimality');
  //let isOptColName = optCol?.name;
  const optColName = '__Optimality';

  form.append(ui.h1('Optimize'));

  numCols.forEach((col, idx) => {
    const feature: NumericFeature = {
      toOptimize: numColsCount - idx - 1 < 2, //false,
      optType: OPT_TYPE.MIN,
    };

    const name = col.name;

    const typeInp = ui.input.choice(name, {
      value: feature.optType,
      nullable: false,
      items: [OPT_TYPE.MIN, OPT_TYPE.MAX],
      onValueChanged: (val) => {
        feature.optType = val;
        computeParetoFront();
        updateViewer();
      },
    });
    ui.tooltip.bind(typeInp.input, 'Type of optimization');

    typeInp.input.hidden = !feature.toOptimize;

    const enableInp = ui.input.toggle('', {
      value: feature.toOptimize,
      onValueChanged: (val) => {
        feature.toOptimize = val;
        typeInp.input.hidden = !val;
        computeParetoFront();
        updateViewer();
      },
    });
    ui.tooltip.bind(enableInp.root, () => `${feature.toOptimize ? 'Dis' : 'En'}able optimization for "${name}"`);

    enableInp.root.classList.add('pareto-switch-input');

    typeInp.root.insertBefore(enableInp.root, typeInp.captionLabel);

    form.append(typeInp.root);
    features.set(name, feature);
  });

  const view = grok.shell.getTableView(df.name);

  view.dockManager.dock(form, DG.DOCK_TYPE.LEFT, null, undefined, 0.25);

  if (numColsCount < 2)
    return;

  const labelCols = cols.filter((col) => {
    if (col.type !== DG.COLUMN_TYPE.STRING)
      return false;

    return col.categories.length === rowCount;
  }).map((col) => col.name);

  const scatter = DG.Viewer.scatterPlot(df, {labelColumnNames: labelCols});
  view.dockManager.dock(scatter, DG.DOCK_TYPE.RIGHT, null, undefined, 0.5);

  const computeParetoFront = () => {
    const data: NumericArray[] = [];
    const sense: OPT_TYPE[] = [];

    features.forEach((fea, name) => {
      if (fea.toOptimize) {
        data.push(df.col(name)!.getRawData());
        sense.push(fea.optType);
      }
    });

    if (data.length > 0) {
      const mask = paretoMaskFromCoordinates(data, sense, rowCount);
      const colOpt = DG.Column.fromStrings('__Optimality', mask);
      df.columns.remove('__Optimality', true);
      df.columns.add(colOpt);
    }
  };

  const updateViewer = () => {
    const colNames: string[] = [];

    features.forEach((fea, name) => {
      if (fea.toOptimize)
        colNames.push(name);
    });

    scatter.setOptions({colorColumnName: '__Optimality'});

    if (colNames.length > 0)
      scatter.setOptions({xColumnName: colNames[0]});

    if (colNames.length > 1)
      scatter.setOptions({yColumnName: colNames[1]});
  }; // updateViewer

  computeParetoFront();
  updateViewer();
} // runParetoOptimizer

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/pareto.css';
import {paretoMaskFromCoordinates} from './pareto-computations';
import {NumericFeature, OPT_TYPE, NumericArray, DIFFERENCE, DOCK_RATIO, OPTIMALITY_COL_NAME} from './defs';


export class ParetoOptimizer {
  private df: DG.DataFrame;
  private numCols: DG.Column[];
  private numColsCount: number;
  private rowCount: number;
  private features = new Map<string, NumericFeature>();
  private view: DG.TableView;
  private scatter: DG.ScatterPlotViewer;
  private scatterNode: DG.DockNode;
  private resultColName: string;
  private labelColumnNames: string[];
  private scatter3dDockNode: DG.DockNode | null = null;

  constructor(df: DG.DataFrame) {
    this.df = df;
    const cols = df.columns;
    const colList = cols.toList();
    this.numCols = colList.filter((col) => col.isNumerical);
    this.numColsCount = this.numCols.length;
    this.rowCount = df.rowCount;
    this.view = grok.shell.getTableView(df.name);
    this.labelColumnNames = this.getLabelColNames();
    this.scatter = DG.Viewer.scatterPlot(df, {labelColumnNames: this.labelColumnNames});
    this.scatterNode = this.view.dockManager.dock(
      this.scatter,
      DG.DOCK_TYPE.RIGHT,
      null,
      '2D Pareto front',
      DOCK_RATIO.VIEWER,
    );

    this.resultColName = cols.getUnusedName(OPTIMALITY_COL_NAME);
  }

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

    this.view.dockManager.dock(form, DG.DOCK_TYPE.LEFT, null, undefined, DOCK_RATIO.FORM);
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
    }
  } // computeParetoFront

  private updateScatter(colNames: string[]): void {
    this.scatter.setOptions({colorColumnName: this.resultColName});

    if (colNames.length > 0)
      this.scatter.setOptions({xColumnName: colNames[0]});

    if (colNames.length > 1)
      this.scatter.setOptions({yColumnName: colNames[1]});
  } // updateScatter

  private build3dScatter(colNames: string[]): void {
    try {
      const plot = DG.Viewer.scatterPlot3d(this.df, {
        xColumnName: colNames[0],
        yColumnName: colNames[1],
        zColumnName: colNames[2],
        colorColumnName: this.resultColName,
      // labelColumnName: this.labelColumnNames[0], // <-- this may crash the page
      });

      this.scatter3dDockNode = this.view.dockManager.dock(
        plot,
        DG.DOCK_TYPE.FILL,
        this.scatterNode,
        '3D Pareto front',
      );
    } catch (err) {
      this.scatter3dDockNode = null;
      grok.shell.warning('Cannot');
    }
  } // build3dScatter

  private updateVisualization(): void {
    const colNames: string[] = [];

    this.features.forEach((fea, name) => {
      if (fea.toOptimize)
        colNames.push(name);
    });

    this.updateScatter(colNames);

    if (this.scatter3dDockNode !== null) {
      this.view.dockManager.close(this.scatter3dDockNode);
      this.scatter3dDockNode = null;
    }

    if (colNames.length > 2)
      this.build3dScatter(colNames);
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

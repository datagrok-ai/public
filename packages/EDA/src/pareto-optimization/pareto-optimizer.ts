import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/pareto.css';
import { NumericArray, Sense, paretoMaskFromCoordinates} from './pareto-computations';

export enum OPT_TYPE {
  MIN = 'Minimize',
  MAX = 'Maximize',
};

export type NumericFeature = {
  toOptimize: boolean,
  optType: OPT_TYPE,
};

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
  let optColName = '__Optimality';

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
  
  let scatter = DG.Viewer.scatterPlot(df, {labelColumnNames: labelCols});  
  view.dockManager.dock(scatter, DG.DOCK_TYPE.RIGHT, null, undefined, 0.5);

  const computeParetoFront = () => {
    const data: NumericArray[] = [];
    const sense: Sense[] = [];

    features.forEach((fea, name) => {
      if (fea.toOptimize) {
        data.push(df.col(name)!.getRawData());
        sense.push(fea.optType === OPT_TYPE.MAX ? "max" : "min");
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

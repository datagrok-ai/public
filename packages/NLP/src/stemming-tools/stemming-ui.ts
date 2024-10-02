/* eslint-disable valid-jsdoc */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MetricInfo, getDefaultMetric, getMetricTypesChoicesList,
  DISTANCE_TYPE, DIST_TYPES_ARR} from './metrics';
import {stemCash, ColUseInfo, getEmbeddings} from './stemming-tools';
import {TEXT_SEM_TYPE, TINY, POLAR_FREQ} from './constants';

/** Modify similarity metric. */
export function modifyMetric(df: DG.DataFrame): void {
  const colsData = new Map<string, ColUseInfo>();
  const colsInp = new Map<string, {metricInput: DG.InputBase, weightInput: DG.InputBase}>();
  const inputElements = new Map<string, HTMLDivElement>();
  const dlg = ui.dialog({title: 'Edit distance'});
  const initCheckedCols = [...stemCash.metricDef!.keys()];

  const onColumnsChanged = (columns: DG.Column[]) => {
    const names = columns.map((col) => col.name);

    distInput.root.hidden = (names.length < 2);
    columnsHeader.hidden = (names.length < 1);

    for (const key of colsData.keys()) {
      const val = colsData.get(key);

      if (names.includes(key)) {
          val!.use = true;
          //@ts-ignore
          inputElements.get(key)?.hidden = false;
      } else {
          val!.use = false;
          //@ts-ignore
          inputElements.get(key)?.hidden = true;
      }

      colsData.set(key, val!);
    }
  };

  const colsInput = ui.input.columns('Features', {value: df.columns.byNames(initCheckedCols), table: df, onValueChanged: (value) => onColumnsChanged(value)});
  colsInput.setTooltip('Features used in computing similarity measure.');
  dlg.add(ui.h3('Source'));
  dlg.add(colsInput);

  const distInput = ui.input.choice('Distance', {value: stemCash.aggrDistance, items: DIST_TYPES_ARR,
    onValueChanged: (value) => {stemCash.aggrDistance = value;}});

  distInput.setTooltip('Type of distance between elements with the specified features.');
  dlg.add(distInput);
  distInput.root.hidden = (initCheckedCols.length < 2);

  const columnsHeader = ui.h3('Features');
  dlg.add(columnsHeader);
  columnsHeader.hidden = (initCheckedCols.length < 1);

  // add an appropriate inputs
  for (const col of df.columns) {
    const name = col.name;

    const colData = {
      type: col.type,
      metric: getDefaultMetric(col),
      use: false,
    };

    if (stemCash.metricDef?.has(name)) {
      const val = stemCash.metricDef.get(name);
      colData.use = true;
      colData.metric.type = val!.type;
      colData.metric.weight = val!.weight;
    }

    colsData.set(name, colData);

    const choices = getMetricTypesChoicesList(col);
    const metricInput = ui.input.choice(`${name}:`, {value: colData.metric.type as string, items: choices,
      onValueChanged: (value) => {
        const val = colsData.get(name);
        //@ts-ignore
        val?.metric.type = value;
        colsData.set(name, val!);
      }});
    metricInput.setTooltip(`Type of metric between the '${name}' feature values.`);

    const weightInput = ui.input.float('metric with the weight', {value: colData.metric.weight,
      onValueChanged: (value) => {
        const val = colsData.get(name);
        //@ts-ignore
        val?.metric.weight = value;
        colsData.set(name, val!);
      }});
    weightInput.setTooltip(`Weight coefficient of the '${name}' feature metric.`);

    const inputs = {metricInput: metricInput, weightInput: weightInput};

    colsInp.set(name, inputs);

    const uiElem = ui.divH([inputs.metricInput.root, inputs.weightInput.root]);
    inputElements.set(name, uiElem);

    uiElem.hidden = !colData.use;

    dlg.add(uiElem);
  } // for

  dlg
    .onCancel(() => {})
    .onOK(() => {
      const map = new Map<string, MetricInfo>();

      for (const key of colsData.keys()) {
        const val = colsData.get(key);

        if (val?.use)
          map.set(key, val.metric);
      }

      stemCash.metricDef = map;
    })
    .show({x: 300, y: 300});
} // modifyMetric

/** Run computation of text embeddings */
export function runTextEmdsComputing(): void {
  const df = grok.shell.t;

  if (df === null) {
    grok.shell.warning('Open dataframe with text column');
    return;
  }

  const columns = df.columns;
  const textColNames = [] as string[];

  for (const col of columns) {
    if (col.semType === TEXT_SEM_TYPE)
    textColNames.push(col.name);
  }

  if (textColNames.length === 0) {
    grok.shell.warning(`No columns with the "${TEXT_SEM_TYPE}" quality`);
    return;
  }

  const computingProps = [
    {"name": "texts", "inputType": "Choice", choices: textColNames, "description": "Name of target text column", "nullable": false},
    {"name": "inPlace", "caption": "In-place", "inputType": "Bool", "description": "Defines whether to add results to the current table or provide them in a new view"},
    {"name": "visualize", "caption": "Visualize", "inputType": "Bool", "description": "Defines whether to add a scatteplot with marked embeddings"},
  ].map((p) => DG.Property.fromOptions(p));  

  const computingOptions = {
    texts: textColNames[0],
    inPlace: true,
    visualize: true,
  };

  const computingOptionsForm = ui.input.form(computingOptions, computingProps);

  const umapProps = [
    {"name": "components", "inputType": "Int", "showPlusMinus": true, min: 2, max: 10, "nullable": false,  "description": "The number of components (dimensions) to project the data to"},
    {"name": "epochs", "inputType": "Int", "showPlusMinus": true, min: 20, max: 800, step: 20, "nullable": false, "description": "The number of epochs to optimize embeddings"},
    {"name": "neighbors", "inputType": "Int", "showPlusMinus": true, min: 2, max: 20, "nullable": false, "description": "The number of nearest neighbors to construct the fuzzy manifold"},
    {"name": "minDist", "caption": "min dist", "inputType": "Float", min: 0.001, max: 2, step: 0.001, "nullable": false, "showSlider": true, "description": "The effective minimum distance between embedded points"},
    {"name": "spread", "inputType": "Float", min: 1, max: 5, step: 1, "showSlider": true, "nullable": false, "description": "The effective scale of embedded points"},
  ].map((p) => DG.Property.fromOptions(p));

  const umapOptions = {
    components: 2,
    epochs: 100,
    neighbors: 4,
    minDist: 0.001,
    spread: 1,
  };

  const umapOptionsForm = ui.input.form(umapOptions, umapProps);

  const acc = ui.accordion();
  const umapPane = acc.addPane('UMAP', () => umapOptionsForm);
  ui.tooltip.bind(umapPane.root, 'Edit UMAP settings');

  const dlg = ui.dialog({title: 'Compute embeddings', helpUrl: '/help/explore/text-embeddings'});

  const runComputations = () => {
    dlg.close();

    try {
      const targetCol: DG.Column = df.col(computingOptions.texts)!;
      const embds = getEmbeddings(df, targetCol, umapOptions);

      if (computingOptions.inPlace) {
        for (const col of embds) {
          col.name = columns.getUnusedName(col.name);
          columns.add(col);
        }
      }
      else {
        const embdsTable = DG.DataFrame.fromColumns([targetCol].concat(embds));
        embdsTable.name = `'${targetCol.name}' embeddings`;
        grok.shell.addTableView(embdsTable);
      }

      if (computingOptions.visualize) {
        const rowCount = df.rowCount;

        const markerSize = new Float32Array(rowCount);
        const markerColor = new Float32Array(rowCount);

        const xMean = embds[0].stats.avg;
        const xStd = embds[0].stats.stdev + TINY; // TINY is added to prevent division by zero
        const xRaw = embds[0].getRawData();
        const yMean = embds[1].stats.avg;
        const yStd = embds[1].stats.stdev + TINY;
        const yRaw = embds[1].getRawData();

        let xNorm: number;
        let yNorm: number;
        let radius: number;
        let angle: number;

        // Marker size & color are specified using polar coordinates
        for (let i = 0; i < rowCount; ++i) {
          // get normalized embeddings
          xNorm = (xRaw[i] - xMean) / xStd;
          yNorm = (yRaw[i] - yMean) / yStd;
          
          // compute polar coordinates
          radius = Math.sqrt(xNorm**2 + yNorm**2);
          angle = Math.acos(xNorm / (radius + TINY)) * (yNorm > 0 ? 1 : -1);

          // heuristics
          markerSize[i] = radius;
          markerColor[i] = Math.sin(1.0 / (TINY + Math.log(radius + TINY))) * Math.sin(POLAR_FREQ * angle);
        }

        const sizeCol = DG.Column.fromFloat32Array('embeddings size', markerSize);
        const colorCol = DG.Column.fromFloat32Array('embeddings color', markerColor);

        sizeCol.name = columns.getUnusedName(sizeCol.name);
        colorCol.name = columns.getUnusedName(colorCol.name);

        const v = grok.shell.v as DG.TableView;
        const t = v.table!;
        t.columns.add(sizeCol);
        t.columns.add(colorCol);
        const visibleColNames = t.columns.names().filter((name) => ![sizeCol.name, colorCol.name].includes(name));
        v.grid.columns.setVisible([visibleColNames[0]]);
        v.grid.columns.setVisible(visibleColNames);

        v.addViewer(DG.VIEWER.SCATTER_PLOT, {
          xColumnName: embds[0].name,
          yColumnName: embds[1].name,
          colorColumnName: colorCol.name,
          sizeColumnName: sizeCol.name,
          jitterSize: 6,
        });
      }
    } 
    catch (error) {
      grok.shell.error((error instanceof Error) ? error.message : 'Core issue');
    }
  };

  dlg.addButton('Run', runComputations, undefined, 'Compute text embeddings using UMAP')

  dlg
    .add(computingOptionsForm)
    .add(acc)
    .show();

  grok.shell.v.append(dlg);
}

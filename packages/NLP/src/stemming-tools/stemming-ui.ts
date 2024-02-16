/* eslint-disable valid-jsdoc */

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MetricInfo, getDefaultMetric, getMetricTypesChoicesList,
  DISTANCE_TYPE, DIST_TYPES_ARR} from './metrics';
import {stemCash, ColUseInfo} from './stemming-tools';

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

  const colsInput = ui.columnsInput('Features', df, onColumnsChanged, {checked: initCheckedCols});
  colsInput.setTooltip('Features used in computing similarity measure.');
  dlg.add(ui.h3('Source'));
  dlg.add(colsInput);

  const distInput = ui.choiceInput(
    'Distance',
    stemCash.aggrDistance,
    DIST_TYPES_ARR,
    (dist: DISTANCE_TYPE) => {stemCash.aggrDistance = dist;},
  );

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
    const metricInput = ui.choiceInput(`${name}:`, colData.metric.type as string, choices, (str: string) => {
      const val = colsData.get(name);
      //@ts-ignore
      val?.metric.type = str;
      colsData.set(name, val!);
    });
    metricInput.setTooltip(`Type of metric between the '${name}' feature values.`);

    const weightInput = ui.floatInput('metric with the weight', colData.metric.weight, (w: number) => {
      const val = colsData.get(name);
      //@ts-ignore
      val?.metric.weight = w;
      colsData.set(name, val!);
    });
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

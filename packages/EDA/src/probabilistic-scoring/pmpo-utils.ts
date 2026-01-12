import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/pmpo.css';

import {COLORS, DESCR_TABLE_TITLE, DESCR_TITLE, DescriptorStatistics, DesirabilityProfileProperties,
  FOLDER, P_VAL, PMPO_COMPUTE_FAILED, PmpoParams, SCORES_TITLE, SELECTED_TITLE, STAT_TO_TITLE_MAP, TINY,
  WEIGHT_TITLE} from './pmpo-defs';
import {computeSigmoidParamsFromX0, getCutoffs, normalPdf, sigmoidS, solveNormalIntersection} from './stat-tools';
import {getColorScaleDiv} from '../pareto-optimization/utils';
import {OPT_TYPE} from '../pareto-optimization/defs';

export function getDescriptorStatisticsTable(stats: Map<string, DescriptorStatistics>): DG.DataFrame {
  const descrCount = stats.size;
  const rawArrs = new Map<string, Float64Array>();

  // Create raw data arrays
  STAT_TO_TITLE_MAP.forEach((_, key) => {
    rawArrs.set(key, new Float64Array(descrCount));
  });

  const descrNames = [...stats.keys()];
  const cols = [DG.Column.fromStrings(DESCR_TITLE, descrNames)];

  // Fill stat columns
  descrNames.forEach((descr, idx) => {
    const curStat = stats.get(descr);

    if (curStat != null) {
      STAT_TO_TITLE_MAP.forEach((_, key) => {
        const val = curStat[key as keyof DescriptorStatistics];
        const arr = rawArrs.get(key);
        arr![idx] = val;
      });
    }
  });

  // Create stat columns
  STAT_TO_TITLE_MAP.forEach((title, field) => {
    cols.push(DG.Column.fromFloat64Array(title, rawArrs.get(field)!));
  });

  // Create the resulting table
  const res = DG.DataFrame.fromColumns(cols);
  res.name = DESCR_TABLE_TITLE;

  return res;
} // getDescriptorStatisticsTable

export function getFilteredByPvalue(descrStats: DG.DataFrame, pValThresh: number): string[] {
  const selected: string[] = [];

  const descrCol = descrStats.col(DESCR_TITLE);

  if (descrCol == null)
    throw new Error(`No column "${DESCR_TITLE} in the table with descriptors statistics.`);

  const descr = descrCol.toList();

  const pValCol = descrStats.col(P_VAL);

  if (pValCol == null)
    throw new Error(`No column "${P_VAL} in the table with descriptors statistics.`);

  const pVals = pValCol.getRawData();

  for (let i = 0; i < descrStats.rowCount; ++i) {
    if (pVals[i] < pValThresh)
      selected.push(descr[i]);
  }

  return selected;
} // getFilteredByPvalue

export function addSelectedDescriptorsCol(descrStats: DG.DataFrame, selected: string[]): DG.DataFrame {
  if (selected.length < 1)
    throw new Error('Empty list of selected descriptors.');

  const rowCount = descrStats.rowCount;
  const selArr = new Array<boolean>(rowCount);
  const descrCol = descrStats.col(DESCR_TITLE);

  if (descrCol == null)
    throw new Error(`No column "${DESCR_TITLE} in the table with descriptors statistics.`);

  const descr = descrCol.toList();
  let res = true;
  const colors: Record<string, string> = {};

  for (let i = 0; i < rowCount; ++i) {
    res = selected.includes(descr[i]);
    selArr[i] = res;
    colors[descr[i]] = res ? COLORS.SELECTED : COLORS.SKIPPED;
  }

  descrCol.colors.setCategorical(colors);

  // Added selected column
  descrStats.columns.add(DG.Column.fromList(DG.COLUMN_TYPE.BOOL, SELECTED_TITLE, selArr));

  return descrStats;
} // addSelectedDescriptorsCol

export function getDescrTooltip(): HTMLElement {
  const firstLine = ui.div();
  firstLine.classList.add('eda-pmpo-tooltip-line');
  const selectedBox = ui.div();
  selectedBox.classList.add('eda-pmpo-box');
  selectedBox.style.backgroundColor = COLORS.SELECTED;
  const selectedLabel = ui.span([]);
  selectedLabel.textContent = '- selected';
  firstLine.appendChild(selectedBox);
  firstLine.appendChild(selectedLabel);

  const secondLine = ui.div();
  secondLine.classList.add('eda-pmpo-tooltip-line');
  const nonSelectedBox = ui.div();
  nonSelectedBox.classList.add('eda-pmpo-box');
  nonSelectedBox.style.backgroundColor = COLORS.SKIPPED;
  const nonSelectedLabel = ui.span([]);
  nonSelectedLabel.textContent = '- excluded';

  secondLine.appendChild(nonSelectedBox);
  secondLine.appendChild(nonSelectedLabel);

  return ui.divV([
    ui.h2(DESCR_TITLE),
    ui.divText('Use of descriptors in model construction:'),
    firstLine,
    secondLine,
  ]);
} // getDescrTooltip

export function getScoreTooltip(): HTMLElement {
  return ui.divV([
    ui.h2(SCORES_TITLE),
    ui.divText('Scores computed using the trained probabilistic multi-parameter optimization (pMPO) model.'),
    getColorScaleDiv(OPT_TYPE.MAX, false),
  ]);
} // getScoreTooltip

export function getFilteredByCorrelations(descriptors: DG.ColumnList, selectedByPvalue: string[],
  statistics: Map<string, DescriptorStatistics>, r2Tresh: number): string[] {
  const getCorrelations = () => {
    const triples: [string, string, number][] = [];

    for (let i = 0; i < selectedByPvalue.length; ++i) {
      for (let j = i + 1; j < selectedByPvalue.length; ++j) {
        triples.push([
          selectedByPvalue[i],
          selectedByPvalue[j],
          descriptors.byName(selectedByPvalue[i]).stats.corr(descriptors.byName(selectedByPvalue[j]))**2,
        ]);
      }
    }

    return triples;
  };

  const correlations = getCorrelations().sort((a, b) => b[2] - a[2]);

  const keep = new Set(selectedByPvalue);

  correlations.filter((triple) => triple[2] > r2Tresh).forEach((triple) => {
    const [descr1, descr2, _] = triple;
    const pVal1 = statistics.get(descr1)!.pValue;
    const pVal2 = statistics.get(descr2)!.pValue;
    const tStat1 = statistics.get(descr1)!.tstat;
    const tStat2 = statistics.get(descr2)!.tstat;

    if (pVal1 > pVal2)
      keep.delete(descr1);
    else if (pVal1 < pVal2)
      keep.delete(descr2);
    else { // process the case of p-value = 0 (or too small)
      if (Math.abs(tStat1) > Math.abs(tStat2))
        keep.delete(descr2);
      else
        keep.delete(descr1);
    }
  });

  return [...keep];
} // getFilteredByCorrelations

export function getModelParams(desired: DG.DataFrame, nonDesired: DG.DataFrame,
  selected: string[], qCutoff: number): Map<string, PmpoParams> {
  const params = new Map<string, PmpoParams>();

  let sum = 0;

  // Compute params
  selected.forEach((name) => {
    const desLen = desired.rowCount;
    const nonDesLen = nonDesired.rowCount;

    const desCol = desired.col(name);
    if (desCol == null)
      throw new Error(PMPO_COMPUTE_FAILED + `: no column "${name}" in the desired table.`);

    const nonDesCol = nonDesired.col(name);
    if (nonDesCol == null)
      throw new Error(PMPO_COMPUTE_FAILED + `: no column "${name}" in the non-desired table.`);

    const muDes = desCol.stats.avg;
    const sigmaDes = desCol.stats.stdev * Math.sqrt((desLen - 1) / desLen);

    const muNonDes = nonDesCol.stats.avg;
    const sigmaNonDes = nonDesCol.stats.stdev * Math.sqrt((nonDesLen - 1) / nonDesLen);

    const cutoffs = getCutoffs(muDes, sigmaDes, muNonDes, sigmaNonDes);
    const intersections = solveNormalIntersection(muDes, sigmaDes, muNonDes, sigmaNonDes);

    let x0: number | null = null;

    if (intersections.length > 0) {
      for (const r of intersections) {
        const low = Math.min(muDes, muNonDes);
        const high = Math.max(muDes, muNonDes);

        if ((low - TINY <= r) && (r <= high + TINY)) {
          x0 = r;
          break;
        }
      }

      if (x0 == null)
        x0 = intersections[0];
    } else
      x0 = cutoffs.cutoff;

    const xBound = cutoffs.cutoffNotDesired;
    const sigmoidParams = computeSigmoidParamsFromX0(muDes, sigmaDes, x0, xBound, qCutoff);

    const z = Math.abs(muDes - muNonDes) / (sigmaDes + sigmaNonDes);
    sum += z;

    params.set(name, {
      desAvg: muDes,
      desStd: sigmaDes,
      nonDesAvg: muNonDes,
      nonDesStd: sigmaNonDes,
      cutoff: cutoffs.cutoff,
      cutoffDesired: cutoffs.cutoffDesired,
      cutoffNotDesired: cutoffs.cutoffNotDesired,
      pX0: sigmoidParams.pX0,
      b: sigmoidParams.b,
      c: sigmoidParams.c,
      zScore: z,
      weight: z,
      intersections: intersections,
      x0: x0,
      xBound: xBound,
    });
  });

  // Normalize weights
  params.forEach((param) => {
    param.weight = param.zScore / sum;
  });

  return params;
} // getModelParams

export function getWeightsTable(params: Map<string, PmpoParams>): DG.DataFrame {
  const count = params.size;
  const descriptors = new Array<string>(count);
  const weights = new Float64Array(count);

  let idx = 0;

  params.forEach((param, name) => {
    descriptors[idx] = name;
    weights[idx] = param.weight;
    ++idx;
  });

  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings(DESCR_TITLE, descriptors),
    DG.Column.fromFloat64Array(WEIGHT_TITLE, weights),
  ]);
} // getWeightsTable

export async function loadPmpoParams(file: DG.FileInfo): Promise<Map<string, PmpoParams>> {
  const jsonText = await file.readAsString();
  const parsedObj = JSON.parse(jsonText);

  return new Map(Object.entries(parsedObj.properties));
} // loadPmpoParams

export async function saveModel(params: Map<string, PmpoParams>, modelName: string): Promise<void> {
  let fileName = modelName;
  const nameInput = ui.input.string('File', {
    value: fileName,
    nullable: false,
    onValueChanged: (val) => {
      fileName = val;
      dlg.getButton('Save').disabled = (fileName.length < 1) || (folderName.length < 1);
    },
  });

  let folderName = FOLDER;
  const folderInput = ui.input.string('Folder', {
    value: folderName,
    nullable: false,
    onValueChanged: (val) => {
      folderName = val;
      dlg.getButton('Save').disabled = (fileName.length < 1) || (folderName.length < 1);
    },
  });

  const save = async () => {
    const path = `${folderName}/${fileName}.json`;
    try {
      const jsonString = JSON.stringify(objectToSave(), null, 2);
      await grok.dapi.files.writeAsText(path, jsonString);
      grok.shell.info(`Saved to ${path}`);
    } catch (err) {
      grok.shell.error(`Failed to save: ${err instanceof Error ? err.message : 'the platform issue'}.`);
    }
    dlg.close();
  };

  const objectToSave = () => {
    if (typeInput.value) {
      return {
        'type': 'MPO Desirability Profile',
        'name': nameInput.value,
        'description': descriptionInput.value,
        'properties': getDesirabilityProfileProperties(params),
      };
    }

    return {
      'type': 'Probabilistic MPO Model',
      'name': nameInput.value,
      'description': descriptionInput.value,
      'properties': Object.fromEntries(params),
    };
  };

  const modelNameInput = ui.input.string('Name', {value: modelName, nullable: true});
  const descriptionInput = ui.input.textArea('Description', {value: ' ', nullable: true});
  const typeInput = ui.input.bool('Desirability Profile', {
    value: true,
    tooltipText: 'Save the model as an MPO Desirability Profile. If disabled, the model is saved in the pMPO format.',
  });

  const dlg = ui.dialog({title: 'Save model'})
    .add(ui.h2('Path'))
    .add(folderInput)
    .add(nameInput)
    .add(ui.h2('Model'))
    .add(modelNameInput)
    .add(descriptionInput)
    .add(typeInput)
    .addButton('Save', async () => {
      const exist = await grok.dapi.files.exists(`${folderName}/${fileName}.json`);
      if (!exist)
        await save();
      else {
        // Handle overwrite confirmation
        ui.dialog({title: 'Warning'})
          .add(ui.label('Overwrite existing file?'))
          .onOK(async () => await save())
          .show();
      }
    })
    .show();
} // saveModel

function getDesirabilityProfileProperties(params: Map<string, PmpoParams>) {
  const props: DesirabilityProfileProperties = {};

  params.forEach((param, name) => {
    const range = significantPoints(param);
    const scale = getScale(param, range);
    props[name] = {
      weight: param.weight * scale,
      line: getLine(param),
      min: Math.min(...range),
      max: Math.max(...range),
    };
  });

  return props;
} // getDesirabilityProfileProperties

function getArgsOfGaussFunc(mu: number, sigma: number): number[] {
  return [
    mu - 3 * sigma,
    mu - 2.5 * sigma,
    mu - 2 * sigma,
    mu - 1.5 * sigma,
    mu - sigma,
    mu - 0.5 * sigma,
    mu - 0.25 * sigma,
    mu,
    mu + 0.25 * sigma,
    mu + 0.5 * sigma,
    mu + sigma,
    mu + 1.5 * sigma,
    mu + 2 * sigma,
    mu + 2.5 * sigma,
    mu + 3 * sigma,
  ];
} // getArgsOfGaussFunc

function getScale(param: PmpoParams, range: number[]): number {
  const values = range.map((x) => basicFunction(x, param));

  return Math.max(...values);
}

function basicFunction(x: number, param: PmpoParams): number {
  return normalPdf(x, param.desAvg, param.desStd) * sigmoidS(x, param.x0, param.b, param.c);
}

function getLine(param: PmpoParams): [number, number][] {
  //const range = getArgsOfGaussFunc(param.desAvg, param.desStd);
  const range = significantPoints(param);
  const scale = getScale(param, range);

  return range.map((x) => [x, basicFunction(x, param) / scale]);
}

function significantPoints(param: PmpoParams): number[] {
  const start = param.desAvg - 10 * param.desStd;
  const end = param.desAvg + 10 * param.desStd;
  const steps = 1000;

  let arg = start;
  let func = basicFunction(arg, param);
  let x = 0;
  let y = 0;

  for (let i = 0; i <= steps; i++) {
    x = start + ((end - start) * i) / steps;
    y = basicFunction(x, param);
    if (y > func) {
      arg = x;
      func = y;
    }
  }

  return getArgsOfGaussFunc(arg, param.desStd);
} // significantPoints

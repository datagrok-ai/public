import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/pmpo.css';

import {COLORS, DESCR_TABLE_TITLE, DESCR_TITLE, DescriptorStatistics, FOLDER, P_VAL, PMPO_COMPUTE_FAILED, PmpoParams,
  SELECTED_TITLE, STAT_TO_TITLE_MAP, TINY, WEIGHT_TITLE} from './pmpo-defs';
import {computeSigmoidParamsFromX0, getCutoffs, solveNormalIntersection} from './stat-tools';

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

  return new Map(Object.entries(parsedObj));
} // loadPmpoParams

export async function saveModel(params: Map<string, PmpoParams>, modelName: string): Promise<void> {
  let fileName = modelName;
  const nameInput = ui.input.string('Name', {
    value: fileName,
    nullable: false,
    onValueChanged: () => {
      if (nameInput.value) {
        fileName = nameInput.value;
        dlg.getButton('Save').disabled = (folderInput.value != null);
      } else
        dlg.getButton('Save').disabled = true;
    },
  });

  let folder = FOLDER;
  const folderInput = ui.input.string('Folder', {
    value: folder,
    nullable: false,
    onValueChanged: () => {
      if (nameInput.value) {
        folder = folderInput.value;
        dlg.getButton('Save').disabled = (nameInput.value != null);
      } else
        dlg.getButton('Save').disabled = true;
    },
  });

  const save = async () => {
    const path = `${folder}/${fileName}.json`;
    try {
      const obj = Object.fromEntries(params);
      const json = JSON.stringify(obj, null, 2);
      await grok.dapi.files.writeAsText(path, json);
      grok.shell.info(`Saved to ${path}`);
    } catch (err) {
      grok.shell.error(`Failed to save: ${err instanceof Error ? err.message : 'the platform issue'}.`);
    }
    dlg.close();
  };

  const dlg = ui.dialog({title: 'Save model'})
    .add(folderInput)
    .add(nameInput)
    .addButton('Save', async () => {
      const exist = await grok.dapi.files.exists(`${folder}/${fileName}.json`);
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

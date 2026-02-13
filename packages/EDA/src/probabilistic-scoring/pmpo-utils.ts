// Utility functions for probabilistic scoring (pMPO)
// Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC4716604/

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/pmpo.css';

import {COLORS, DESCR_TABLE_TITLE, DESCR_TITLE, DescriptorStatistics, DesirabilityProfileProperties,
  DESIRABILITY_COL_NAME, FOLDER, P_VAL, PMPO_COMPUTE_FAILED, PmpoParams, SCORES_TITLE,
  SELECTED_TITLE, STAT_TO_TITLE_MAP, TINY, WEIGHT_TITLE, CorrelationTriple,
  BASIC_RANGE_SIGMA_COEFFS, EXTENDED_RANGE_SIGMA_COEFFS} from './pmpo-defs';
import {computeSigmoidParamsFromX0, getCutoffs, gaussDesirabilityFunc, sigmoidS,
  solveNormalIntersection} from './stat-tools';
import {getColorScaleDiv} from '../pareto-optimization/utils';
import {OPT_TYPE} from '../pareto-optimization/defs';

/** Returns a DataFrame with descriptor statistics.
 * @param stats Map of descriptor names to their statistics.
 */
export function getDescriptorStatisticsTable(stats: Map<string, DescriptorStatistics>): DG.DataFrame {
  const descrCount = stats.size;
  const rawArrs = new Map<string, Float64Array>();

  // Create raw data arrays
  STAT_TO_TITLE_MAP.forEach((_, key) => {
    rawArrs.set(key, new Float64Array(descrCount));
  });

  const descrNames = [...stats.keys()];
  const cols = [
    DG.Column.fromStrings(DESCR_TITLE, descrNames),
    DG.Column.fromInt32Array(DESIRABILITY_COL_NAME, new Int32Array(descrCount)),
    DG.Column.fromFloat32Array(WEIGHT_TITLE, new Float32Array(descrCount).fill(DG.FLOAT_NULL)),
  ];

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

/** Returns names of descriptors with p-value below the given threshold.
 * @param descrStats DataFrame with descriptor statistics.
 * @param pValThresh P-value threshold.
 */
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

/** Adds a boolean column indicating whether each descriptor is selected.
 * @param descrStats DataFrame with descriptor statistics.
 * @param selected List of selected descriptor names.
 */
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

/** Returns tooltip element describing descriptor selection colors. */
export function getDescrTooltip(title: string, text: string, selected: string, excluded: string): HTMLElement {
  const firstLine = ui.div();
  firstLine.classList.add('eda-pmpo-tooltip-line');
  const selectedBox = ui.div();
  selectedBox.classList.add('eda-pmpo-box');
  selectedBox.style.backgroundColor = COLORS.SELECTED;
  const selectedLabel = ui.span([]);
  selectedLabel.textContent = `- ${selected}`;
  firstLine.appendChild(selectedBox);
  firstLine.appendChild(selectedLabel);

  const secondLine = ui.div();
  secondLine.classList.add('eda-pmpo-tooltip-line');
  const nonSelectedBox = ui.div();
  nonSelectedBox.classList.add('eda-pmpo-box');
  nonSelectedBox.style.backgroundColor = COLORS.SKIPPED;
  const nonSelectedLabel = ui.span([]);
  nonSelectedLabel.textContent = `- ${excluded}`;

  secondLine.appendChild(nonSelectedBox);
  secondLine.appendChild(nonSelectedLabel);

  return ui.divV([ui.h2(title), ui.divText(text), firstLine, secondLine]);
} // getDescrTooltip

/** Returns tooltip element describing score colors. */
export function getScoreTooltip(): HTMLElement {
  return ui.divV([
    ui.h2(SCORES_TITLE),
    ui.divText('Scores computed using the trained probabilistic multi-parameter optimization (pMPO) model.'),
    getColorScaleDiv(OPT_TYPE.MAX, false),
  ]);
} // getScoreTooltip

/** Returns list of descriptor correlation triples.
 * @param descriptors Descriptor column list.
 * @param selectedByPvalue List of descriptor names selected by p-value.
 */
export function getCorrelationTriples(descriptors: DG.ColumnList, selectedByPvalue: string[]): CorrelationTriple[] {
  const triples: CorrelationTriple[] = [];

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
} // getCorrelationTriples

/** Returns names of descriptors filtered by correlation threshold.
 * @param descriptors Descriptor column list.
 * @param selectedByPvalue List of descriptor names selected by p-value.
 */
export function getFilteredByCorrelations(descriptors: DG.ColumnList, selectedByPvalue: string[],
  statistics: Map<string, DescriptorStatistics>, r2Tresh: number, correlationTriples: CorrelationTriple[]): string[] {
  const correlations = correlationTriples.sort((a, b) => b[2] - a[2]);

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

/** Computes pMPO model parameters for selected descriptors.
 * @param desired DataFrame with desired compounds.
 * @param nonDesired DataFrame with non-desired compounds.
 * @param selected List of selected descriptor names.
 * @param qCutoff Q-value cutoff.
 */
export function getModelParams(desired: DG.DataFrame, nonDesired: DG.DataFrame,
  selected: string[], qCutoff: number): Map<string, PmpoParams> {
  const params = new Map<string, PmpoParams>();

  let sum = 0;

  // Compute params for each selected descriptor
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

    // Unbiased standard deviation
    const sigmaDes = desCol.stats.stdev * Math.sqrt((desLen - 1) / desLen);

    const muNonDes = nonDesCol.stats.avg;

    // Unbiased standard deviation
    const sigmaNonDes = nonDesCol.stats.stdev * Math.sqrt((nonDesLen - 1) / nonDesLen);

    // Compute cutoffs
    const cutoffs = getCutoffs(muDes, sigmaDes, muNonDes, sigmaNonDes);

    // column_stats['inflection'] = np.exp(-np.square((column_stats['cutoff'] - column_stats['good_mean'])) /
    //                                     (2 * np.square(column_stats['good_std'])))
    // Compute inflection point
    const inflection = Math.exp(-((cutoffs.cutoff - muDes) ** 2) / (2 * (sigmaDes ** 2)));

    // Compute intersections of the two normal distributions
    const intersections = solveNormalIntersection(muDes, sigmaDes, muNonDes, sigmaNonDes);

    const b = (Math.pow(inflection, -1.0) - 1.0);
    const n = (Math.pow(qCutoff, -1.0) - 1.0);
    const c = Math.pow(10.0, ((Math.log10(n / b)) / (-1.0 * (muNonDes - cutoffs.cutoff))));

    // Compute parameters for the generalized sigmoid function TODO: delete

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

    // Store computed parameters
    params.set(name, {
      desAvg: muDes,
      desStd: sigmaDes,
      nonDesAvg: muNonDes,
      nonDesStd: sigmaNonDes,
      min: Math.min(desCol.stats.min, nonDesCol.stats.min),
      max: Math.max(desCol.stats.max, nonDesCol.stats.max),
      cutoff: cutoffs.cutoff,
      cutoffDesired: cutoffs.cutoffDesired,
      cutoffNotDesired: cutoffs.cutoffNotDesired,
      pX0: sigmoidParams.pX0,
      b: b,
      c: c,
      zScore: z,
      weight: z,
      intersections: intersections,
      x0: x0,
      xBound: xBound,
      inflection: inflection,
    });
  });

  // Normalize weights
  params.forEach((param) => {
    param.weight = param.zScore / sum;
  });

  return params;
} // getModelParams

/** Returns a DataFrame with descriptor weights.
 * @param params Map of descriptor names to their pMPO parameters.
 */
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

/** Loads pMPO model parameters from a file.
 * @param file FileInfo object pointing to the JSON model file.
 */
export async function loadPmpoParams(file: DG.FileInfo): Promise<Map<string, PmpoParams>> {
  const jsonText = await file.readAsString();
  const parsedObj = JSON.parse(jsonText);

  return new Map(Object.entries(parsedObj.properties));
} // loadPmpoParams

/** Returns JSON object representing an MPO Desirability Profile.
 * @param params Map of descriptor names to their pMPO parameters.
 * @param name Name of the desirability profile.
 * @param description Description of the desirability profile.
 */
export function getDesirabilityProfileJson(params: Map<string, PmpoParams>, useSigmoidalCorrection: boolean,
  name: string, description: string, truncatedRange: boolean): any {
  return {
    'type': 'MPO Desirability Profile',
    'name': name,
    'description': description,
    'properties': getDesirabilityProfileProperties(params, useSigmoidalCorrection, truncatedRange),
  };
}

/** Saves pMPO model parameters to a file.
 * @param params Map of descriptor names to their pMPO parameters.
 * @param modelName Suggested model name (used as default file name).
 */
export async function saveModel(params: Map<string, PmpoParams>, modelName: string,
  useSigmoidalCorrection: boolean): Promise<void> {
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
      return getDesirabilityProfileJson(
        params,
        useSigmoidalCorrection,
        nameInput.value,
        descriptionInput.value,
        false,
      );
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

/** Adds columns with correlation coefficients between descriptors.
 * @param df DataFrame to which the columns will be added.
 * @param descriptorNames List of descriptor names.
 * @param triples List of descriptor correlation triples.
 * @param selectedByCorr List of descriptor names selected after correlation filtering.
 */
export function addCorrelationColumns(df: DG.DataFrame, descriptorNames: string[],
  triples: CorrelationTriple[], selectedByCorr: string[]): DG.DataFrame {
  const raw = new Map<string, Float32Array>();
  descriptorNames.forEach((name) => {
    raw.set(name, new Float32Array(descriptorNames.length).fill(DG.FLOAT_NULL));
  });

  const descrColVals = df.col(DESCR_TITLE)!.toList();

  triples.forEach((triple) => {
    const [descr1, descr2, r2] = triple;

    raw.get(descr1)![descrColVals.indexOf(descr2)] = r2;
    raw.get(descr2)![descrColVals.indexOf(descr1)] = r2;
    raw.get(descr1)![descrColVals.indexOf(descr1)] = 1;
    raw.get(descr2)![descrColVals.indexOf(descr2)] = 1;
  });

  selectedByCorr.forEach((name) => {
    df.columns.add(DG.Column.fromFloat32Array(name, raw.get(name)!));
  });

  return df;
} // addCorrelationColumns

/** Sets color coding for the p-value column in the statistics table
 * @param table DataFrame with descriptor statistics.
 * @param pValTresh P-value threshold.
*/
export function setPvalColumnColorCoding(table: DG.DataFrame, pValTresh: number): void {
  const pValCol = table.col(P_VAL);
  if (pValCol == null)
    return;

  const rules: Record<string, string> = {};
  rules[`<${pValTresh}`] = COLORS.SELECTED;
  rules[`>=${pValTresh}`] = COLORS.SKIPPED;

  pValCol.meta.colors.setConditional(rules);
} // setPvalColumnColorCoding

/** Sets color coding for the correlation columns in the statistics table.
 * @param table DataFrame with descriptor statistics.
 * @param descriptorNames List of descriptor names.
 * @param r2Tresh R-squared threshold.
*/
export function setCorrColumnColorCoding(table: DG.DataFrame, descriptorNames: string[], r2Tresh: number): void {
  descriptorNames.forEach((name) => {
    const col = table.col(name);
    if (col == null)
      return;

    const rules: Record<string, string> = {};
    rules[`>=${r2Tresh}`] = COLORS.SKIPPED;
    rules[`<${r2Tresh}`] = COLORS.SELECTED;

    col.meta.colors.setConditional(rules);
  });
} // setCorrColumnColorCoding

/** Returns desirability profile properties for the given pMPO parameters.
 * @param params Map of descriptor names to their pMPO parameters.
 * @param useSigmoidalCorrection Whether to use sigmoidal correction in desirability functions.
 * @param displayProfile Whether to create a profile to be displayed in the stat grid (true - truncated range).
 */
function getDesirabilityProfileProperties(params: Map<string, PmpoParams>,
  useSigmoidalCorrection: boolean, truncatedRange: boolean): DesirabilityProfileProperties {
  const props: DesirabilityProfileProperties = {};

  params.forEach((param, name) => {
    const range = significantPoints(param, truncatedRange);
    props[name] = {
      weight: param.weight,
      line: getLine(param, useSigmoidalCorrection, truncatedRange),
      min: Math.min(...range),
      max: Math.max(...range),
    };
  });

  return props;
} // getDesirabilityProfileProperties

/** Returns array of arguments for Gaussian function centered at mu with stddev sigma.
 * @param mu Mean of the Gaussian function.
 * @param sigma Standard deviation of the Gaussian function.
 * @param truncatedRange Whether to use truncated range (for interactive app) or extended range (for full profile).
 * @return Array of arguments for the Gaussian function.
*/
function getArgsOfGaussFunc(mu: number, sigma: number, truncatedRange: boolean): number[] {
  return truncatedRange ?
    BASIC_RANGE_SIGMA_COEFFS.map((coeff) => mu + coeff * sigma) : // range for interactive app
    EXTENDED_RANGE_SIGMA_COEFFS.map((coeff) => mu + coeff * sigma); // actual full range for desirability profile
} // getArgsOfGaussFunc

/** Basic pMPO function combining Gaussian and sigmoid functions.
 * @param x Argument.
 * @param param pMPO parameters.
 * @param useSigmoidalCorrection Whether to use sigmoidal correction.
 * @return Value of the basic pMPO function at x.
*/
function basicFunction(x: number, param: PmpoParams, useSigmoidalCorrection: boolean): number {
  return gaussDesirabilityFunc(x, param.desAvg, param.desStd) *
    (useSigmoidalCorrection ? sigmoidS(x, param.cutoff, param.b, param.c) : 1);
}

/** Returns line points for the given pMPO parameters.
 * @param param pMPO parameters.
 * @param useSigmoidalCorrection Whether to use sigmoidal correction.
 * @param truncatedRange Whether to use truncated range (for interactive app) or extended range (for full profile).
 * @return Array of [x, y] points representing the desirability function line.
*/
function getLine(param: PmpoParams, useSigmoidalCorrection: boolean, truncatedRange: boolean): [number, number][] {
  const range = significantPoints(param, truncatedRange);

  return range.map((x) => [x, basicFunction(x, param, useSigmoidalCorrection)]);
}

/** Returns significant points for the given pMPO parameters.
 * @param param pMPO parameters.
 * @param truncatedRange Whether to use truncated range (for interactive app) or extended range (for full profile).
 * @return Array of significant points for the desirability function.
*/
function significantPoints(param: PmpoParams, truncatedRange: boolean): number[] {
  const points = getArgsOfGaussFunc(param.desAvg, param.desStd, truncatedRange);

  /* Truncate range to show less points */
  if (truncatedRange) {
    const min = Math.min(param.min, param.desAvg - 3 * param.desStd);
    const max = Math.max(param.max, param.desAvg + 3 * param.desStd);

    return points
      .filter((x) => (min <= x) && (x <= max))
      .sort();
  }

  return points;
} // significantPoints

/** Custom error class for pMPO-related errors. */
export class PmpoError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'PmpoError';
  }
}

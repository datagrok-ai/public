/* eslint-disable camelcase */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from './package';
import type {BuilderType} from '../escher_src/src/Builder';
import type {CobraModelData, FluxHistogram, ReactionBounds, SamplingFunctionResult} from '../escher_src/src/ts/types';
import {WorkerCobraSolver} from './cobra';
import type {PrecomputedWarmup} from './cobra/sampler-wrapper';
import {ItemsGrid} from '@datagrok-libraries/utils/src/items-grid';
import {solveUsingGLPKJvail} from './cobra/glpkJS';
import {openTimeCourseDialog, clearTimeCourseSlider} from './timeCourse';

export function parsePath(path?: string): string | null {
  if (!path)
    return null;
  path = path.replaceAll('%20', ' ');
  if (path.startsWith('/'))
    return path.substring(1);
  return path;
}

export type MetabolicAnalysisState = {
  description?: string,
} & Record<string, any>;

/** Name of the analysis last saved/loaded per builder — the single source the UI dialogs,
 * the URL path loader, and the AI view functions all report and update. */
const currentAnalysisByBuilder = new WeakMap<BuilderType, string>();

export function getCurrentAnalysis(builder: BuilderType): string | null {
  return currentAnalysisByBuilder.get(builder) ?? null;
}

export async function saveStateDialog(builder: BuilderType, currentCampaign?: string) {
  currentCampaign ??= getCurrentAnalysis(builder) ?? undefined;
  const dialog = ui.dialog('Save Analysis');
  const existingCampaigns: string[] = [];
  let currentExists = false;
  if (currentCampaign && await _package.files.exists(`campaigns/${currentCampaign}.json`))
    currentExists = true;

  const saveTypeInput = ui.input.choice('Save To',
    {items: ['New Analysis', 'Existing Analysis'], value: currentExists ? 'Existing Analysis' : 'New Analysis'});
  const nameInput = ui.input.string('Name', {nullable: true, placeholder: 'Analysis name'});
  const descriptionInput = ui.input.textArea('Description', {nullable: true});
  nameInput.addValidator((v) => existingCampaigns.includes(v) ? 'Analysis with this name already exists' : null);
  const campList = await loadCampaigns();
  if (campList)
    existingCampaigns.push(...campList);

  const campaignChoicInput = ui.input.choice('Name', {items: existingCampaigns, nullable: false});
  if (existingCampaigns.length < 1) {
    saveTypeInput.value = 'New Analysis';
    saveTypeInput.enabled = false;
    ui.tooltip.bind(saveTypeInput.root, 'No existing analyses found');
    campaignChoicInput.root.style.display = 'none';
    nameInput.root.style.display = 'flex';
  }

  saveTypeInput.onChanged.subscribe((_) => {
    if (saveTypeInput.value === 'New Analysis') {
      campaignChoicInput.root.style.display = 'none';
      nameInput.root.style.display = 'flex';
    } else {
      campaignChoicInput.root.style.display = 'flex';
      nameInput.root.style.display = 'none';
    }
  });

  if (currentExists) {
    saveTypeInput.value = 'Existing Analysis';
    campaignChoicInput.value = currentCampaign!;
  }

  campaignChoicInput.onChanged.subscribe(async (_) => {
    try {
      if (campaignChoicInput.value) {
        const camp = existingCampaigns.find((c) => c === campaignChoicInput.value);
        if (camp) {
          const campJSON = JSON.parse(await _package.files.readAsText(`campaigns/${camp}.json`));
          descriptionInput.value = campJSON.description ?? '';
        }
      }
    } catch (e) {
      grok.shell.error('Error loading analysis');
      console.error(e);
    }
  });

  saveTypeInput.fireChanged();

  dialog.add(saveTypeInput);
  dialog.add(nameInput);
  dialog.add(campaignChoicInput);
  dialog.add(descriptionInput);
  dialog.addButton('Save', async () => {
    dialog.close();
    const name = saveTypeInput.value === 'New Analysis' ? nameInput.value : campaignChoicInput.value;
    if (!name) {
      grok.shell.warning('Name is required for analysis');
      return;
    }
    await saveAnalysisState(builder, name, descriptionInput.value);
    grok.shell.info(`Analysis ${name} saved`);
  });
  dialog.show();
  dialog.root.addEventListener('keydown', (e) => {
    //needed so that keystrokes are not recognised by escher
    e.stopPropagation();
  });
}

function replaceUrlPath(curAnalysisName: string) {
  const curHref = window.location.href.toLowerCase();
  const lastIndex = curHref.lastIndexOf('metabolicgraph');
  if (lastIndex < 0)
    return;
  const newHref = curHref.substring(0, lastIndex) + `metabolicgraph/${curAnalysisName}`;
  if (curHref !== newHref) {
    // @ts-ignore
    if (history.replaceState) {
      const title = document.title;
      const obj = {Title: title, Url: newHref};
      history.replaceState(obj, obj.Title, obj.Url);
    }
  }
}

export async function loadCampaigns() {
  // ui.setUpdateIndicator(dialog.root, true);
  const pg = DG.TaskBarProgressIndicator.create('Loading analysis list');
  try {
    const camps = (await _package.files.list('campaigns')).filter((f) => f.name.endsWith('.json'))
      .map((f) => f.name.substring(0, f.name.length - 5));
    pg.close();
    return camps;
  } catch (e) {
    grok.shell.error('Error loading analysis list');
    console.error(e);
  }
  pg.close();
  return null;
}

export function analysisExists(name: string): Promise<boolean> {
  return _package.files.exists(`campaigns/${name}.json`);
}

export function readAnalysis(name: string): Promise<string> {
  return _package.files.readAsText(`campaigns/${name}.json`);
}

export async function saveAnalysisState(builder: BuilderType, name: string, description?: string | null) {
  const state = builder.getSavingState();
  state.description = description ?? '';
  await _package.files.writeAsText(`campaigns/${name}.json`, JSON.stringify(state));
  currentAnalysisByBuilder.set(builder, name);
  replaceUrlPath(name);
}

export async function loadAnalisisDialog(builder: BuilderType) {
  const campList = await loadCampaigns();
  if (!campList || campList.length < 1) {
    grok.shell.warning('No analyses found');
    return;
  }
  const dialog = ui.dialog('Load Analysis');
  const campaignChoicInput = ui.input.choice('Name', {items: campList, nullable: false});
  dialog.add(campaignChoicInput);
  dialog.addButton('Load', async () => {
    dialog.close();
    const camp = campList.find((c) => c === campaignChoicInput.value);
    if (camp)
      loadStateProxy(builder, await readAnalysis(camp), camp);
  });
  dialog.show();
}

export async function loadStateProxy(builder: BuilderType, campJSON: string, path?: string) {
  builder.loadSavingState(campJSON);
  if (path)
    currentAnalysisByBuilder.set(builder, path);
  replaceUrlPath(path ?? '');
  builder.settings.set('loadAction', () => loadAnalisisDialog(builder));
  builder.settings.set('saveAction', () => saveStateDialog(builder, path));
}

/** cobra's own OptGP warmup points, computed on the server (ComputeExtremePoints.py). */
export async function computeWarmupPython(cobraModel: CobraModelData): Promise<PrecomputedWarmup> {
  const func = DG.Func.find({package: 'MetabolicGraph', name: 'ComputeExtremePoints'})[0];
  if (!func)
    throw new Error('ComputeExtremePoints Python function not found');
  const modelString = JSON.stringify(cobraModel).replaceAll(`'{}'`, '{}').replaceAll('"{}"', '{}');
  const resultString: string = await func.apply({cobraModel: modelString});
  return JSON.parse(resultString) as PrecomputedWarmup;
}

export async function sampleReactions(cobraModel: CobraModelData, builder: BuilderType): Promise<SamplingFunctionResult> {
  return new Promise<SamplingFunctionResult>((resolve, reject) => {
    const samplesInput = ui.input.int('Number of samples', {value: 10000, nullable: false});
    const binsInput = ui.input.int('Number of bins', {value: 20, nullable: false, tooltipText: 'Number of bins for histogram in the reaction tooltips'});
    const thinningInput = ui.input.int('Thinning', {value: 10, nullable: false, tooltipText: 'Thinning interval for the sampler'});
    const aggregationInput = ui.input.choice<FluxAggregation>('Aggregation', {items: FLUX_AGGREGATIONS, value: 'Mean',
      nullable: false, tooltipText: 'How each reaction\'s samples are reduced to the single flux shown on the map. ' +
        'Mean keeps the fluxes mass-balanced, Median resists skewed distributions, Mode is the histogram peak'});
    const addDfInput = ui.input.bool('Add DataFrame', {value: false, tooltipText: 'Add a DataFrame with all sampled fluxes to the workspace'});
    const runUsingPythonInput = ui.input.bool('Run using Python', {value: false, tooltipText: 'Use Python OptGpSampling script instead of WebAssembly sampler'});
    const usePythonFBAInput = ui.input.bool('Use Python FBA', {value: false, tooltipText: 'Compute the warmup points (cobra\'s minimum and maximum of every flux) with cobra in Python, then sample using WebAssembly. The browser computes the same points itself, so this only adds a server round trip'});
    const getInput = () => {
      return {samples: samplesInput.value, thinning: thinningInput.value, bins: binsInput.value,
        aggregation: aggregationInput.value ?? 'Mean', addDf: addDfInput.value, runInPython: runUsingPythonInput.value,
        usePythonFBA: usePythonFBAInput.value};
    };
    const applyInput = (x: Partial<ReturnType<typeof getInput>>) => {
      x.samples && (samplesInput.value = x.samples);
      x.thinning && (thinningInput.value = x.thinning);
      x.bins && (binsInput.value = x.bins);
      x.aggregation && FLUX_AGGREGATIONS.includes(x.aggregation) && (aggregationInput.value = x.aggregation);
      x.addDf != undefined && (addDfInput.value = x.addDf);
      x.runInPython != undefined && (runUsingPythonInput.value = x.runInPython);
      x.usePythonFBA != undefined && (usePythonFBAInput.value = x.usePythonFBA);
    };
    const innerLocalStorageKey = 'metabolic-graph-reaction-sampling-dialog-inner-last-input';
    runUsingPythonInput.onChanged.subscribe(() => {
      if (runUsingPythonInput.value) {
        usePythonFBAInput.value = false;
        usePythonFBAInput.enabled = false;
      } else
        usePythonFBAInput.enabled = true;
    });

    const dlg = ui.dialog('Sample Reactions')
      .add(samplesInput)
      .add(thinningInput)
      .add(binsInput)
      .add(aggregationInput)
      .add(addDfInput)
      .add(runUsingPythonInput)
      .add(usePythonFBAInput);
    dlg.addButton('Time-course…', () => {
      // hand off to the interpolated-bounds flow, reusing the parameters set above;
      // resolve the single-reaction sampling promise as cancelled so its distribution is left untouched
      resolve({data: new Map(), cancled: true});
      dlg.close();
      openTimeCourseDialog(cobraModel, builder, getInput());
    });
    dlg.onOK(async () => {
      try {
        try {
          const toSave = getInput();
          localStorage.setItem(innerLocalStorageKey, JSON.stringify(toSave));
        } catch (_) {}
        let precomputedWarmup: PrecomputedWarmup | undefined;
        if (usePythonFBAInput.value) {
          const pg = DG.TaskBarProgressIndicator.create('Computing warmup points with cobra (Python)');
          try {
            precomputedWarmup = await computeWarmupPython(cobraModel);
          } finally {
            pg.close();
          }
        }
        const aggregation = aggregationInput.value ?? 'Mean';
        const result = runUsingPythonInput.value ?
          await runReactionSamplingPython(cobraModel, builder, binsInput.value, addDfInput.value, samplesInput.value, thinningInput.value, aggregation) :
          await runReactionSampling(cobraModel, builder, binsInput.value, addDfInput.value, samplesInput.value, thinningInput.value, aggregation, precomputedWarmup);
        resolve(result);
      } catch (e) {
        grok.shell.error('Error sampling reactions');
        console.error(e);
        reject(e);
      }
    })
      .onCancel(() => {
        resolve({data: new Map(), cancled: true});
      })
      .show({center: true})
      .history(() => getInput(),
        (x) => {
          applyInput(x);
        });
    try {
      const lastInput = localStorage.getItem(innerLocalStorageKey);
      if (lastInput) {
        const parsed = JSON.parse(lastInput);
        applyInput(parsed);
      }
    } catch (_) {}
  });
}

/**
 * Statistic that reduces a reaction's samples to the one flux shown on the map: the central statistics
 * of the COBRA Toolbox's calcSampleStats. Harmonic/geometric means are left out: fluxes are signed and
 * often zero.
 */
export type FluxAggregation = 'Mean' | 'Median' | 'Mode';
export const FLUX_AGGREGATIONS: FluxAggregation[] = ['Mean', 'Median', 'Mode'];

/** 'median flux' etc., for data source captions. */
export const aggregatedFluxLabel = (aggregation: FluxAggregation) => `${aggregation.toLowerCase()} flux`;

/** With every reaction fixed at zero there is nothing to sample. */
const allFluxesFixedAtZero = (cobraModel: CobraModelData) =>
  cobraModel.reactions.every((r) => (r.lower_bound ?? 0) >= 0 && (r.upper_bound ?? 0) <= 0);

/**
 * Bin one reaction's samples over their own [min, max] and reduce them to a single flux.
 * Takes ownership of `values`: the median sorts it in place.
 */
function summarizeFluxes(values: Float32Array | Float64Array, bins: number, aggregation: FluxAggregation):
  {histogram: FluxHistogram, flux: number} | null {
  const n = values.length;
  if (n === 0)
    return null;
  let min = Infinity;
  let max = -Infinity;
  let sum = 0;
  for (let i = 0; i < n; i++) {
    const v = values[i];
    if (v < min) min = v;
    if (v > max) max = v;
    sum += v;
  }
  // a spread within float noise is one fixed flux: a single bin at its value
  if (max - min <= 1e-6 * Math.max(1, Math.abs(min), Math.abs(max))) {
    min = max = (min + max) / 2;
    bins = 1;
  }
  const counts = new Array<number>(bins).fill(0);
  const binWidth = (max - min) / bins;
  if (bins === 1)
    counts[0] = n;
  else {
    for (let i = 0; i < n; i++)
      counts[Math.min(bins - 1, Math.floor((values[i] - min) / binWidth))]++;
  }

  let flux = sum / n;
  if (aggregation === 'Median') {
    values.sort();
    const mid = n >> 1;
    flux = n % 2 ? values[mid] : (values[mid - 1] + values[mid]) / 2;
  } else if (aggregation === 'Mode') {
    // center of the fullest bin: the peak of the tooltip histogram
    let top = 0;
    for (let b = 1; b < bins; b++) {
      if (counts[b] > counts[top])
        top = b;
    }
    flux = min + (top + 0.5) * binWidth;
  }
  return {histogram: {min, max, counts}, flux};
}

/** Per-reaction histograms (tooltips) and aggregated fluxes (map colors); `column(i)` returns a fresh copy. */
function reduceSamples(ids: string[], column: (i: number) => Float32Array | Float64Array, bins: number,
  aggregation: FluxAggregation) {
  const reactionData: {[key: string]: number} = {};
  const distribution: SamplingFunctionResult = {data: new Map()};
  const binCount = Math.max(1, Math.floor(bins));
  ids.forEach((id, i) => {
    const summary = summarizeFluxes(column(i), binCount, aggregation);
    if (summary) {
      reactionData[id] = summary.flux;
      distribution.data.set(id, summary.histogram);
    }
  });
  return {reactionData, distribution};
}

/** Build a DataFrame (one column per reaction) from raw stacked flux samples. */
export function samplesToDataFrame(cobraModel: CobraModelData, results: Float32Array, nSamples: number): DG.DataFrame {
  const reactions = cobraModel.reactions;
  const columns = reactions.map((r, j) => {
    const col = DG.Column.float(r.id, nSamples);
    col.init((i) => results[i * reactions.length + j]);
    return col;
  });
  return DG.DataFrame.fromColumns(columns);
}

export type FluxSamplingResult = {
  reactionData: {[key: string]: number}; // per-reaction aggregated flux (drives map colors)
  distribution: SamplingFunctionResult; // per-reaction histograms over their own ranges (drives tooltips)
  results?: Float32Array; // raw samples (WASM path), stacked row-wise
  df?: DG.DataFrame; // raw samples as a DataFrame (Python path)
};

/** Run the WebAssembly sampler and reduce to aggregated fluxes + histograms. No UI/builder side effects. */
export async function sampleFluxesWasm(cobraModel: CobraModelData, bins = 20, nSamples = 1000, thinning = 20,
  aggregation: FluxAggregation = 'Mean', precomputedWarmup?: PrecomputedWarmup): Promise<FluxSamplingResult | null> {
  if (allFluxesFixedAtZero(cobraModel))
    return null;
  const results = await WorkerCobraSolver.runSampling(cobraModel, nSamples, thinning, precomputedWarmup);
  const reactionCount = cobraModel.reactions.length;
  const {reactionData, distribution} = reduceSamples(cobraModel.reactions.map((r) => r.id), (i) => {
    const values = new Float32Array(nSamples);
    for (let j = 0; j < nSamples; j++)
      values[j] = results[j * reactionCount + i];
    return values;
  }, bins, aggregation);
  return {reactionData, distribution, results};
}

/** Run the Python OptGpSampling function and reduce to aggregated fluxes + histograms. No UI/builder side effects. */
export async function sampleFluxesPython(cobraModel: CobraModelData, bins = 20, nSamples = 1000, thinning = 10,
  aggregation: FluxAggregation = 'Mean'): Promise<FluxSamplingResult | null> {
  const func = DG.Func.find({package: 'MetabolicGraph', name: 'OptGpSampling'})[0];
  if (!func) {
    grok.shell.error('OptGpSampling function not found');
    return null;
  }
  if (allFluxesFixedAtZero(cobraModel))
    return null;
  const modelString = JSON.stringify(cobraModel).replaceAll(`'{}'`, '{}').replaceAll('"{}"', '{}');
  const resDf: DG.DataFrame = await func.apply({cobraModel: modelString, nSamples: nSamples, thinning: thinning});
  const {reactionData, distribution} = reduceSamples(resDf.columns.names(), (i) => {
    const col = resDf.columns.byIndex(i);
    // the raw buffer can outlive the column length, and marks missing values in-band
    const values = Float64Array.from(col.getRawData().subarray(0, col.length));
    return col.stats.missingValueCount ? values.filter((_, r) => !col.isNone(r)) : values;
  }, bins, aggregation);
  return {reactionData, distribution, df: resDf};
}

export async function runReactionSampling(cobraModel: CobraModelData, builder: BuilderType, bins: number = 20,
  addDataFrame: boolean = false, nSamples: number = 1000, thinning: number = 20, aggregation: FluxAggregation = 'Mean',
  precomputedWarmup?: PrecomputedWarmup): Promise<SamplingFunctionResult> {
  if (!cobraModel)
    throw new Error('Cannot run optimization without a model loaded');
  clearTimeCourseSlider(); // a fresh single-shot sampling supersedes any time-course animation
  const pg = DG.TaskBarProgressIndicator.create('Sampling reactions');
  try {
    const res = await sampleFluxesWasm(cobraModel, bins, nSamples, thinning, aggregation, precomputedWarmup);
    if (!res) {
      grok.shell.error('Invalid reaction bounds in the model. Check reaction lower and upper bounds');
      return {data: new Map()};
    }
    builder.set_reaction_data(res.reactionData, `Sampling, ${aggregatedFluxLabel(aggregation)}`);
    if (addDataFrame) {
      const table = samplesToDataFrame(cobraModel, res.results!, nSamples);
      table.name = 'Sampler results';
      grok.shell.addTableView(table);
    }
    return res.distribution;
  } catch (e) {
    grok.shell.error('Error sampling reactions');
    console.error(e);
    return {data: new Map()};
  } finally {
    pg.close();
  }
}

export async function runReactionSamplingPython(cobraModel: CobraModelData, builder: BuilderType, bins: number = 20,
  addDataFrame: boolean = false, nSamples: number = 1000, thinning: number = 10,
  aggregation: FluxAggregation = 'Mean'): Promise<SamplingFunctionResult> {
  clearTimeCourseSlider(); // a fresh single-shot sampling supersedes any time-course animation
  const pg = DG.TaskBarProgressIndicator.create('Sampling reactions (Python)');
  try {
    const res = await sampleFluxesPython(cobraModel, bins, nSamples, thinning, aggregation);
    if (!res) {
      grok.shell.error('Invalid reaction bounds in the model. Check reaction lower and upper bounds');
      return {data: new Map()};
    }
    builder.set_reaction_data(res.reactionData, `Sampling (Python), ${aggregatedFluxLabel(aggregation)}`);
    if (addDataFrame && res.df) {
      res.df.name = 'Sampler results (Python)';
      grok.shell.addTableView(res.df);
    }
    return res.distribution;
  } catch (e) {
    grok.shell.error('Error sampling reactions');
    console.error(e);
    return {data: new Map()};
  } finally {
    pg.close();
  }
}

export function handleReactionDataUpload(view: DG.ViewBase, builder: BuilderType) {
  // // @ts-ignore
  // window.currentBuilder = builder; // for debugging purposes, to access the builder from the console

  // // @ts-ignore
  // window.solver = WorkerCobraSolver;

  view.subs.push(grok.events.onFileImportRequest.subscribe(async (ff) => {
    const f = ff as unknown as DG.EventData<DG.FileImportArgs>;
    // the current view can be the JS view itself or the Dart host wrapping it
    const cur = grok.shell.v as any;
    const isCurrent = cur != null && (cur === view || cur.jsView === view);
    if (!isCurrent || (!f.args.file?.name?.endsWith('.json') && !f.args.file?.name?.endsWith('.csv')))
      return;
    f.preventDefault();
    clearTimeCourseSlider(); // dropping new data supersedes any time-course animation
    const fileString = await f.args.file.text();
    // check JSON
    if (f.args.file.name?.endsWith('.json')) {
      try {
        const parsed = JSON.parse(fileString);
        if (parsed && Object.values(parsed).every((v) => typeof v === 'number')) {
          // this is a valid reaction data JSON
          builder.set_reaction_data(parsed, f.args.file.name);
          grok.shell.info(`Reaction flux data from ${f.args.file.name} loaded`);
        } else
          throw new Error('Invalid reaction data JSON');
      } catch (e) {
        const df = DG.DataFrame.fromJson(fileString);
        grok.shell.addTableView(df);
      }
    } else {
      const df = DG.DataFrame.fromCsv(fileString);
      try {
        // if it is a csv, it can be one of two things: reaction flux data or reaction bounds data. the second one has 3 columns
        const isFluxData = df.columns.length === 2 && df.columns.byIndex(1).isNumerical && df.columns.byIndex(0).isCategorical;
        const isBoundsData = df.columns.length === 3 && df.columns.byIndex(1).isNumerical && df.columns.byIndex(2).isNumerical &&
          df.columns.byIndex(0).isCategorical;
        if (isFluxData) {
          const reactionData: {[key: string]: number} = {};
          const reactionCol = df.columns.byIndex(0);
          const fluxCol = df.columns.byIndex(1);
          for (let i = 0; i < df.rowCount; i++) {
            const reaction = reactionCol.get(i) as string;
            const flux = fluxCol.get(i) as number;
            if (typeof reaction !== 'string' || typeof flux !== 'number')
              continue;

            reactionData[reaction] = flux;
          }
          builder.set_reaction_data(reactionData, f.args.file.name);
          grok.shell.info(`Reaction flux data from ${f.args.file.name} loaded`);
        } else if (isBoundsData) {
          const reactionBoundsData: {[key: string]: ReactionBounds} = {};
          const reactionCol = df.columns.byIndex(0);
          const lowerBoundCol = df.columns.byIndex(1);
          const upperBoundCol = df.columns.byIndex(2);
          for (let i = 0; i < df.rowCount; i++) {
            const reaction = reactionCol.get(i) as string;
            const lowerBound = lowerBoundCol.get(i) as number;
            const upperBound = upperBoundCol.get(i) as number;
            if (typeof reaction !== 'string' || typeof lowerBound !== 'number' || typeof upperBound !== 'number')
              continue;

            reactionBoundsData[reaction] = {lower_bound: lowerBound, upper_bound: upperBound};
          }
          builder.set_reaction_bounds(reactionBoundsData);
          grok.shell.info(`Reaction bounds data from ${f.args.file.name} loaded`);
        } else
          throw new Error('Invalid reaction data CSV');
      } catch (e) {
        grok.shell.addTableView(df);
        return;
      }
    }
  }));
}

export async function runFBADialog(builder: BuilderType) {
  const dialog = ui.dialog('Run FBA');
  const reactions = builder.model_data.reactions;
  const objectiveReactions = reactions.filter((r) => !!r.objective_coefficient);
  const reactionProp = DG.Property.fromOptions({name: 'Reaction', type: DG.TYPE.STRING, choices: reactions.map((r) => r.id), nullable: true});
  const objectiveProp = DG.Property.fromOptions({name: 'Objective', type: DG.TYPE.STRING, choices: ['Maximize', 'Minimize'], nullable: false, defaultValue: 'Maximize'});

  const existingItems = objectiveReactions.map((r) => ({
    Reaction: r.id,
    Objective: r.objective_coefficient ? (r.objective_coefficient > 0 ? 'Maximize' : 'Minimize') : 'Maximize',
  }));

  const propGrid = new ItemsGrid([reactionProp, objectiveProp], existingItems, {addButtonTooltip: 'Add reaction objective', removeButtonTooltip: 'Remove reaction objective', newItemFunction: () => ({Reaction: null, Objective: 'Maximize'})});
  dialog.add(propGrid.root);

  dialog.addButton('Run', async () => {
    clearTimeCourseSlider(); // running FBA supersedes any time-course animation
    const validItems = [...propGrid.items, propGrid.addingItem].filter((i) => i.Reaction && i.Objective);
    // go through existing reactions and delete all the objective coefficients
    for (const r of reactions) {
      if (r.objective_coefficient)
        delete r.objective_coefficient;
      const p = validItems.find((i) => i.Reaction === r.id);
      if (p) {
        // set the objective coefficient
        r.objective_coefficient = p.Objective === 'Maximize' ? 1 : -1;
      }
    }

    dialog.close();
    // run FBA
    try {
      console.log(await solveUsingGLPKJvail(builder.model_data));
      // const solution = await runGLPKFBA(builder.model_data);
      const solution = await WorkerCobraSolver.run_optimization(builder.model_data);
      if (!solution) {
        grok.shell.error('FBA failed to run');
        return;
      }
      const fluxes: Record<string, number> = {};
      solution.reactionNames.forEach((r, i) => {
        fluxes[r] = solution.fluxes[i];
      });
      builder.set_reaction_data(fluxes, 'FBA Result: ' + validItems.map((i) => `${i.Objective} ${i.Reaction}`).join(', '));
    } catch (e) {
      grok.shell.error('FBA failed to run');
      console.error(e);
      return;
    }
  });

  dialog.show();
  dialog.root.style.minWidth = '400px';
}

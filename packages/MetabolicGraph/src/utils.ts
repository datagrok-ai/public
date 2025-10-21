/* eslint-disable camelcase */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from './package';
import type {BuilderType} from '../escher_src/src/Builder';
import type {CobraModelData, ReactionBounds, SamplingFunctionResult} from '../escher_src/src/ts/types';
import {WorkerCobraSolver} from './cobra';
import {ItemsGrid} from '@datagrok-libraries/utils/src/items-grid';
import {solveUsingGLPKJvail} from './cobra/glpkJS';

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

export async function saveStateDialog(builder: BuilderType, currentCampaign?: string) {
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
    const state = builder.getSavingState();
    state.description = descriptionInput.value;
    await _package.files.writeAsText(`campaigns/${name}.json`, JSON.stringify(state));
    grok.shell.info(`Analysis ${name} saved`);
    replaceUrlPath(name);
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

async function loadCampaigns() {
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
    if (camp) {
      const campJSON = await _package.files.readAsText(`campaigns/${camp}.json`);
      loadStateProxy(builder, campJSON, camp);
    }
  });
  dialog.show();
}

export async function loadStateProxy(builder: BuilderType, campJSON: string, path?: string) {
  builder.loadSavingState(campJSON);
  replaceUrlPath(path ?? '');
  builder.settings.set('loadAction', () => loadAnalisisDialog(builder));
  builder.settings.set('saveAction', () => saveStateDialog(builder, path));
}

export async function sampleReactions(cobraModel: CobraModelData, builder: BuilderType): Promise<SamplingFunctionResult> {
  return new Promise<SamplingFunctionResult>((resolve, reject) => {
    const samplesInput = ui.input.int('Number of samples', {value: 10000, nullable: false});
    const binsInput = ui.input.int('Number of bins', {value: 20, nullable: false, tooltipText: 'Number of bins for histogram in the reaction tooltips'});
    const thinningInput = ui.input.int('Thinning', {value: 10, nullable: false, tooltipText: 'Thinning interval for the sampler'});
    const addDfInput = ui.input.bool('Add DataFrame', {value: false, tooltipText: 'Add a DataFrame with all sampled fluxes to the workspace'});
    const runUsingPythonInput = ui.input.bool('Run using Python', {value: false, tooltipText: 'Use Python OptGpSampling script instead of WebAssembly sampler'});
    const getInput = () => {
      return {samples: samplesInput.value, thinning: thinningInput.value, bins: binsInput.value, addDf: addDfInput.value, runInPython: runUsingPythonInput.value};
    };
    const applyInput = (x: Partial<ReturnType<typeof getInput>>) => {
      x.samples && (samplesInput.value = x.samples);
      x.thinning && (thinningInput.value = x.thinning);
      x.bins && (binsInput.value = x.bins);
      x.addDf != undefined && (addDfInput.value = x.addDf);
      x.runInPython != undefined && (runUsingPythonInput.value = x.runInPython);
    };
    const innerLocalStorageKey = 'metabolic-graph-reaction-sampling-dialog-inner-last-input';
    ui.dialog('Sample Reactions')
      .add(samplesInput)
      .add(thinningInput)
      .add(binsInput)
      .add(addDfInput)
      .add(runUsingPythonInput)
      .onOK(async () => {
        try {
          try {
            const toSave = getInput();
            localStorage.setItem(innerLocalStorageKey, JSON.stringify(toSave));
          } catch (_) {}
          const result = runUsingPythonInput.value ?
            await runReactionSamplingPython(cobraModel, builder, samplesInput.value, thinningInput.value) :
            await runReactionSampling(cobraModel, builder, binsInput.value, addDfInput.value, samplesInput.value, thinningInput.value);
          resolve(result);
        } catch (e) {
          grok.shell.error('Error sampling reactions');
          console.error(e);
          reject(e);
        }
      })
      .onCancel(() => {
        resolve({upper_bound: 10, lower_bound: -10, data: new Map<string, number[]>(), cancled: true});
      })
      .show()
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

export async function runReactionSampling(cobraModel: CobraModelData, builder: BuilderType, bins: number = 20, addDataFrame: boolean = false, nSamples: number = 1000, thinning: number = 20): Promise<SamplingFunctionResult> {
  if (!cobraModel)
    throw new Error('Cannot run optimization without a model loaded');

  const pg = DG.TaskBarProgressIndicator.create('Sampling reactions');
  const sampleMap = new Map<string, number[]>();
  const lower_bound = cobraModel.reactions.reduce((min, r) => Math.min(min, r.lower_bound ?? 0), 0);
  const upper_bound = cobraModel.reactions.reduce((max, r) => Math.max(max, r.upper_bound ?? 0), 0);
  if (lower_bound >= upper_bound) {
    grok.shell.error('Invalid reaction bounds in the model. Check reaction lower and upper bounds');
    pg.close();
    return {upper_bound: 10, lower_bound: -10, data: sampleMap};
  }

  try {
    // first run FBA to get a feasible solution
    const results = await WorkerCobraSolver.runSampling(cobraModel, nSamples, thinning);
    const range = upper_bound - lower_bound;
    const step = range / bins;
    const reactionCount = cobraModel.reactions.length;

    for (let i = 0; i < reactionCount; i++)
      sampleMap.set(cobraModel.reactions[i].id, new Array(bins).fill(0));


    const reactionAverages = new Float32Array(reactionCount).fill(0);
    const reactionData: {[key: string]: number} = {};

    for (let i = 0; i < cobraModel.reactions.length; i++) {
      for (let j = 0; j < nSamples; j++) {
        const indexInResult = j * reactionCount + i;
        const flux = results[indexInResult];
        const binIndex = Math.min(bins - 1, Math.max(0, Math.floor((flux - lower_bound) / step)));
        const reactionId = cobraModel.reactions[i].id;
        sampleMap.get(reactionId)![binIndex]++;
        reactionAverages[i] += flux;
      }
    }
    // calculate averages
    for (let i = 0; i < cobraModel.reactions.length; i++) {
      const reactionId = cobraModel.reactions[i].id;
      reactionAverages[i] /= nSamples;
      reactionData[reactionId] = reactionAverages[i];
    }

    builder.set_reaction_data(reactionData, 'Sampling Histogram');

    if (addDataFrame) {
      const reactions = cobraModel.reactions;
      const columns = reactions.map((r, j) => {
        const col = DG.Column.float(r.id, nSamples);
        col.init((i) => results[i * reactions.length + j]);
        return col;
      });
      const table = DG.DataFrame.fromColumns(columns);
      table.name = 'Sampler results';
      grok.shell.addTableView(table);
    }
  } catch (e) {
    grok.shell.error('Error sampling reactions');
    console.error(e);
  }


  pg.close();
  return {upper_bound, lower_bound, data: sampleMap};
}

/// Deprecated. Use runReactionSampling instead
export async function runReactionSamplingPython(cobraModel: CobraModelData, builder: BuilderType, nSamples: number = 1000, thinning: number = 10): Promise<SamplingFunctionResult> {
  const func = DG.Func.find({package: 'MetabolicGraph', name: 'OptGpSampling'})[0];
  const sampleMap = new Map<string, number[]>();
  let lower_bound = -10;
  let upper_bound = 10;
  if (!func) {
    grok.shell.error('OptGpSampling function not found');
    return {upper_bound, lower_bound, data: sampleMap};
  }
  const pg = DG.TaskBarProgressIndicator.create('Sampling reactions');
  try {
    const modelString = JSON.stringify(cobraModel).replaceAll(`'{}'`, '{}').replaceAll('"{}"', '{}');
    const resDf: DG.DataFrame = await func.apply({cobraModel: modelString, nSamples: nSamples, thinning: thinning});
    for (const col of resDf.columns) {
      if (col.name?.toLowerCase() === 'bin range')
        continue;
      sampleMap.set(col.name, col.toList().slice(0, resDf.rowCount - 1)); // get only the first nSamples values
    }
    // there is a column with the bounds called 'Bin Range'
    // parse the column and get the bounds. its format is '[-10, 10)'
    const boundsCol = resDf.col('Bin Range');
    if (boundsCol) {
      const lower: string = boundsCol.get(0);
      const upper: string = boundsCol.get(boundsCol.length - 2); // -2 because the last value is the average, so we need the second last
      lower_bound = Math.floor(parseFloat(lower.trim().substring(1, lower.indexOf(','))));
      upper_bound = Math.ceil(parseFloat(upper.trim().substring(upper.indexOf(',') + 1, upper.length - 1)));
    }

    // average results using samples for each reaction to get average flux
    const reactionData: {[key: string]: number} = {};
    for (const col of resDf.columns) {
      if (col.name?.toLowerCase() === 'bin range')
        continue;
      const reaction = col.name;
      // the last row is the average, so we need to put that into the reactionData
      const averageFlux = col.get(resDf.rowCount - 1);
      reactionData[reaction] = averageFlux;
    }
    builder.set_reaction_data(reactionData, 'Samling Average');
  } catch (e) {
    grok.shell.error('Error sampling reactions');
    console.error(e);
  }
  // calculate average reaction flux for each reaction
  pg.close();
  return {upper_bound, lower_bound, data: sampleMap};
}

export function handleReactionDataUpload(view: DG.View, builder: BuilderType) {
  // // @ts-ignore
  // window.currentBuilder = builder; // for debugging purposes, to access the builder from the console

  // // @ts-ignore
  // window.solver = WorkerCobraSolver;

  grok.events.onFileImportRequest.subscribe(async (ff) => {
    const f = ff as unknown as DG.EventData<DG.FileImportArgs>;
    if ((grok.shell.v as DG.View)?.id !== view.id || (!f.args.file?.name?.endsWith('.json') && !f.args.file?.name?.endsWith('.csv')))
      return;
    f.preventDefault();
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
  });
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

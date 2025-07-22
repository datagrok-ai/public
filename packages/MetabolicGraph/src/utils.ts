/* eslint-disable camelcase */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from './package';
import type {BuilderType} from '../escher_src/src/Builder';
import type {CobraModelData} from '../escher_src/src/ts/types';

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

export async function sampleReactions(cobraModel: CobraModelData, builder: BuilderType) {
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
    const nSamples = 1000;
    const modelString = JSON.stringify(cobraModel).replaceAll(`'{}'`, '{}').replaceAll('"{}"', '{}');
    const resDf: DG.DataFrame = await func.apply({cobraModel: modelString, nSamples: nSamples});
    for (const col of resDf.columns)
      sampleMap.set(col.name, col.toList());
    // there is a column with the bounds called 'Bin Range'
    // parse the column and get the bounds. its format is '[-10, 10)'
    const boundsCol = resDf.col('Bin Range');
    if (boundsCol) {
      const lower: string = boundsCol.get(0);
      const upper: string = boundsCol.get(boundsCol.length - 1);
      lower_bound = Math.floor(parseFloat(lower.trim().substring(1, lower.indexOf(','))));
      upper_bound = Math.ceil(parseFloat(upper.trim().substring(upper.indexOf(',') + 1, upper.length - 1)));
    }

    // average results using samples for each reaction to get average flux
    const bins = resDf.rowCount;
    const range = upper_bound - lower_bound;
    const singleBinRange = range / bins;
    const binMidPoints = Array.from({length: bins}, (_, i) => lower_bound + (i + 0.5) * singleBinRange);
    const reactionData: {[key: string]: number} = {};
    for (const col of resDf.columns) {
      if (col.name?.toLowerCase() === 'bin range')
        continue;
      const reaction = col.name;
      const values = col.toList();
      let sum = 0;
      for (let i = 0; i < bins; i++)
        sum += values[i] * binMidPoints[i];
      reactionData[reaction] = sum / nSamples;
    }
    builder.set_reaction_data(reactionData);
  } catch (e) {
    grok.shell.error('Error sampling reactions');
    console.error(e);
  }
  // calculate average reaction flux for each reaction
  pg.close();
  return {upper_bound, lower_bound, data: sampleMap};
}

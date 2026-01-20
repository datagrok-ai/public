/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {_package} from '../../package';
import * as api from '../../package-api';
import {AppName, HitDesignCampaign} from '../types';
import {CampaignTableName, HitDesignMolColName, ViDColName} from '../consts';


export function obfuscateSmiles(smiles: string): string {
  const SECRET_KEY = _package.meta?.sok;
  if (!SECRET_KEY)
    throw new Error('Secret key for SMILES obfuscation is not set');
  let result = '';
  for (let i = 0; i < smiles.length; i++) {
    const char = smiles.charCodeAt(i);
    const keyChar = SECRET_KEY.charCodeAt(i % SECRET_KEY.length);
    result += String.fromCharCode(char ^ keyChar);
  }
  return btoa(result); // Browser's base64 encode
}

export function deobfuscateSmiles(obfuscated: string): string {
  const SECRET_KEY = _package.meta?.sok;
  if (!SECRET_KEY)
    throw new Error('Secret key for SMILES deobfuscation is not set');
  const decoded = atob(obfuscated); // Browser's base64 decode
  let result = '';
  for (let i = 0; i < decoded.length; i++) {
    const char = decoded.charCodeAt(i);
    const keyChar = SECRET_KEY.charCodeAt(i % SECRET_KEY.length);
    result += String.fromCharCode(char ^ keyChar);
  }
  return result;
}

function getCurrentUserName() {
  const user = DG.User.current();
  return user?.login ?? user?.friendlyName ?? user.name ?? 'unknown user';
}

/**
   * Adds a molecule to the dictionary (or not if it exists) and returns its V-iD string
   * @param mol - molecule string (assumed that it is not canonical at this point)
   */
export async function registerMol(mol: string, campaignID: string, appName: AppName): Promise<string> {
  const canonical = _package.convertToSmiles(mol);
  if (!canonical || canonical === '' || canonical === 'MALFORMED_INPUT_VALUE')
    return '';
  const obfuscated = obfuscateSmiles(canonical);
  const result = await api.queries.addMolecule(obfuscated, appName, campaignID, getCurrentUserName());
  if (result.rowCount === 0 || !result.col('vid')) {
    grok.shell.error('Molecule registration returned no data');
    return '';
  }
  return result.col('vid')!.getString(0);
}

export async function registerMolsBatch(mols: string[], campaignID: string, appName: AppName): Promise<Map<string, string>> {
  const molToVidMap = new Map<string, string>(); // the key should be original mol string, not the canonical one
  // const molSet = new Set<string>(mols.filter((m) => !!m));
  // const molList = Array.from(molSet);

  const canonicalMolsMap = new Map<string, number[]>(); // maps canonical smiles to list of indexes in original array
  mols.forEach((m, i) => {
    const can = _package.convertToSmiles(m);
    if (!can || can === '' || can === 'MALFORMED_INPUT_VALUE')
      return;
    if (!canonicalMolsMap.has(can))
      canonicalMolsMap.set(can, []);
    canonicalMolsMap.get(can)!.push(i); // store original indexes, not just one because we might have duplicates where its same molecule but one is smiles, and other is smth else :D
  });

  if (canonicalMolsMap.size === 0)
    return molToVidMap;

  const canonicalMolList = Array.from(canonicalMolsMap.keys());
  const batchSize = 50;
  for (let i = 0; i < canonicalMolList.length; i += batchSize) {
    const batch = canonicalMolList.slice(i, i + batchSize);
    if (batch.length === 0)
      continue;
    try {
      const obfuscatedBatch = batch.map((m) => obfuscateSmiles(m));
      const result = await api.queries.addMolecules(obfuscatedBatch, appName, campaignID, getCurrentUserName());
      if (!result || result.rowCount !== batch.length || !result.col('vid')) {
        console.error('Molecule registration returned incomplete data for batch', `${i}-${i + batch.length}`, `Batch size: ${batch.length}, Result rows: ${result?.rowCount}. Vid column found: ${!!result?.col('vid')}`);
        throw new Error('Resulting molecule registration has no or incomplete data');
      }
      // at this point we are SURE that all molecules in the batch are unique
      const vidCol = result.col('vid')!.toList();
      for (let j = 0; j < result.rowCount; j++) {
        const originalIndexes = canonicalMolsMap.get(batch[j])!; // original index
        const vid = vidCol[j];
        for (const originalIndex of originalIndexes)
          molToVidMap.set(mols[originalIndex], vid);
      }
    } catch (e) {
      grok.shell.error('Error registering molecules batch. see console for details.');
      _package.logger.error(e);
    }
  }
  return molToVidMap;
}

/**
 * This function will go through all campaign molecules and register all of them into vid storage
 * It is like a migration function - needs to be done ONLY ONCE
 */
export async function registerAllCampaignMols() {
  const appName: AppName = 'Hit Design';
  const allCampaigns = await _package.loadCampaigns(appName);
  const campaignNames = Object.keys(allCampaigns);
  const campaignsSelector = ui.input.choice('Select Campaign', {value: null, nullable: true, items: campaignNames, tooltipText: 'Select campaign to register its molecules. Leave empty to register all campaigns'});

  const toProcess = await new Promise<HitDesignCampaign[] | null>(async (resolve) => {
    ui.dialog('Register Existing Campaign Molecules')
      .add(ui.divText('This operation will go through all molecules in the selected campaign(s) and register them into V-iD storage.\n This operation is irreversible.'))
      .add(ui.divText('You can select a specific campaign or leave it empty to process all campaigns.'))
      .add(ui.divText('After this opteration, all VIDs throughout all campaigns will become unique, and old V-iDs will no longer be valid.'))
      .add(ui.divText('WARNING: Depending on the number of molecules, this operation can take a long time.'))
      .add(campaignsSelector)
      .onOK(() => {
        if (campaignsSelector.value) {
          const selectedCampaign = allCampaigns[campaignsSelector.value];
          resolve([selectedCampaign]);
        } else {
          const allSelectedCampaigns = campaignNames.map((name) => allCampaigns[name]);
          resolve(allSelectedCampaigns);
        }
      }).onCancel(() => resolve(null))
      .show();
  });

  const pg = DG.TaskBarProgressIndicator.create('Registering molecules...');
  if (!toProcess || toProcess.length === 0) {
    pg.close();
    return;
  }

  const total = toProcess.length;
  const count = 0;
  for (const campaign of toProcess) {
    pg.update((count / total)* 100, `Processing campaign: ${campaign.name}`);
    try {
      const fileLoc = `System.AppData/HitTriage/${appName}/campaigns`;
      const campaignTablePath = campaign?.savePath ?? `${fileLoc}/${campaign.name}/${CampaignTableName}`; // p.s. campaign.name is its unique id
      const table = await grok.dapi.files.readCsv(campaignTablePath);
      if (!table) { // might be old campaign, so just skip
        grok.shell.warning(`Failed to load campaign table for campaign ${campaign.name} at path ${campaignTablePath}`);
        console.warn(`Failed to load campaign table for campaign ${campaign.name} at path ${campaignTablePath}`);
        continue;
      }
      const molCol = table.col(HitDesignMolColName);
      if (!molCol) {
        grok.shell.warning(`Molecule column not found in campaign ${campaign.name} table`);
        console.warn(`Molecule column not found in campaign ${campaign.name} table`);
        continue;
      }
      const vidCol = table.col(ViDColName);
      if (!vidCol) {
        grok.shell.warning(`V-iD column not found in campaign ${campaign.name} table`);
        console.warn(`V-iD column not found in campaign ${campaign.name} table`);
        continue;
      }

      const mols = molCol.toList();
      const molsToVid = await registerMolsBatch(mols as string[], campaign.name, appName);
      for (let i = 0; i < table.rowCount; i++) {
        const mol = mols[i];
        if (!mol || mol === '') {
          vidCol.set(i, '', false);
          continue;
        }
        const vid = molsToVid.get(mol);
        if (vid)
          vidCol.set(i, vid, false); // commit only on last set
      }
      vidCol.fireValuesChanged();
      // finally, save the updated table back
      await grok.dapi.files.writeAsText(campaignTablePath, table.toCsv());
      grok.shell.info(`Finished processing campaign ${campaign.name}`);
      console.log(`Finished processing campaign ${campaign.name}`);
    } catch (e) {
      grok.shell.error(`Error processing campaign ${campaign.name}. See console for details.`);
      console.error(`Error processing campaign ${campaign.name} molecules registration:`, e);
    }
  }
  pg.close();
}

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {IKnimeClient} from './knime-client';
import {CachedDeploymentEntry} from './types';
import {registerFuncFromSpec, buildParamMeta, buildSignature, sanitizeFuncName} from './function-registry';

const STORAGE_NAME = 'KnimeLinkFuncCache';
const META_KEY = '_meta';

export function loadCachedEntries(): CachedDeploymentEntry[] {
  const all = grok.userSettings.get(STORAGE_NAME);
  if (!all || !all[META_KEY])
    return [];

  let ids: string[];
  try { ids = JSON.parse(all[META_KEY]); }
  catch { return []; }

  const entries: CachedDeploymentEntry[] = [];
  for (const id of ids) {
    const raw = all[id];
    if (!raw)
      continue;
    try { entries.push(JSON.parse(raw)); }
    catch { /* skip corrupted entry */ }
  }
  return entries;
}

function saveCacheEntries(entries: CachedDeploymentEntry[]): void {
  const data: {[key: string]: string} = {};
  data[META_KEY] = JSON.stringify(entries.map((e) => e.deployment.id));
  for (const entry of entries)
    data[entry.deployment.id] = JSON.stringify(entry);
  grok.userSettings.put(STORAGE_NAME, data);
}

export function registerFromCache(entries: CachedDeploymentEntry[], client: IKnimeClient): void {
  for (const entry of entries) {
    try {
      registerFuncFromSpec(entry.deployment, entry.spec, client);
    }
    catch (e) {
      console.warn(`KnimeLink: failed to register cached function for ${entry.deployment.name}:`, e);
    }
  }
}

export async function refreshAndUpdateCache(client: IKnimeClient): Promise<void> {
  const deployments = await client.listDeployments('rest');

  const oldAll = grok.userSettings.get(STORAGE_NAME) ?? {};
  const oldEntries = new Map<string, CachedDeploymentEntry>();
  try {
    const ids: string[] = oldAll[META_KEY] ? JSON.parse(oldAll[META_KEY]) : [];
    for (const id of ids) {
      if (oldAll[id]) {
        try { oldEntries.set(id, JSON.parse(oldAll[id])); }
        catch { /* skip */ }
      }
    }
  }
  catch { /* no previous cache */ }

  const freshEntries: CachedDeploymentEntry[] = [];

  for (const dep of deployments) {
    try {
      const spec = await client.getWorkflowInputs(dep.id);
      const funcName = sanitizeFuncName(dep.name, dep.id);
      const paramMeta = buildParamMeta(spec.inputs);
      const signature = buildSignature(funcName, paramMeta, spec.outputs);

      const old = oldEntries.get(dep.id);
      if (!old || old.signature !== signature)
        registerFuncFromSpec(dep, spec, client);

      freshEntries.push({deployment: dep, spec, signature});
    }
    catch (e) {
      console.warn(`KnimeLink: failed to fetch spec for ${dep.name}, keeping cached version:`, e);
      const old = oldEntries.get(dep.id);
      if (old)
        freshEntries.push(old);
    }
  }

  saveCacheEntries(freshEntries);
}

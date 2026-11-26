import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Column, ENGINE, ListedModel, ModelReference, NamedModel, OwnRow, Registration, SyncDiff, TYPES}
  from './types';

function columnsOf(json?: string): Column[] {
  try {
    const columns: unknown = JSON.parse(json ?? '');
    return Array.isArray(columns) ?
      columns.map((c: {name: string; type: string; semType?: string}) =>
        ({...c, type: TYPES[c.type] ?? DG.COLUMN_TYPE.STRING})) : [];
  } catch {
    return [];
  }
}

function targetOf(model: ListedModel): Column {
  return columnsOf(model.outputs)[0] ?? {name: 'prediction', type: DG.COLUMN_TYPE.FLOAT};
}

function isModel(listed: ListedModel): listed is NamedModel {
  return !!listed.name;
}

function isFailedConnection(listed: ListedModel): boolean {
  return !isModel(listed) && !!listed.error;
}

function referenceOf(row: OwnRow): Partial<ModelReference> {
  return JSON.parse(row.Options);
}

export function syncDiff(listed: ListedModel[], rows: OwnRow[]): SyncDiff {
  const keyOf = (connectionId?: string, name?: string) => `${connectionId}/${name}`;
  const byKey = new Map<string, OwnRow>();
  const remove: OwnRow[] = [];
  for (const row of rows) {
    const {connectionId, name} = referenceOf(row);
    const key = keyOf(connectionId, name);
    if (byKey.has(key))
      remove.push(row);
    else
      byKey.set(key, row);
  }
  const register: Registration[] = [];
  for (const model of listed.filter(isModel)) {
    const key = keyOf(model.connection_id, model.name);
    const row = byKey.get(key);
    byKey.delete(key);
    if (model.error)
      continue;
    if (row === undefined || isChanged(model, row))
      register.push({model, row});
  }
  const failed = new Set(listed.filter(isFailedConnection).map((m) => m.connection_id));
  for (const row of byKey.values()) {
    if (!failed.has(referenceOf(row).connectionId ?? ''))
      remove.push(row);
  }
  return {register, remove};
}

// A registered version never changes, so the feature names are enough to see a new signature.
function isChanged(model: NamedModel, row: OwnRow): boolean {
  const target = targetOf(model);
  const {version, type} = referenceOf(row);
  return version !== (model.version ?? '') || type !== target.type || row.Target !== target.name ||
    (row.Features ?? '') !== columnsOf(model.inputs).map((c) => c.name).join(', ');
}

async function listedModels(): Promise<ListedModel[]> {
  return (await grok.functions.call('mlflowModels', {})).toJson();
}

async function ownRows(): Promise<OwnRow[]> {
  return (await grok.functions.call('Forge:listModels', {method: ENGINE, own: true})).toJson();
}

function registerModel({model, row}: Registration): Promise<string> {
  const target = targetOf(model);
  const reference: ModelReference = {connectionId: model.connection_id, name: model.name,
    version: model.version ?? '', type: target.type};
  const isClass = target.type === DG.COLUMN_TYPE.STRING || target.type === DG.COLUMN_TYPE.BOOL;
  // The name and the description of a registered row are the user's.
  return grok.functions.call('Forge:registerModel', {id: row?.Id, name: row?.Name ?? model.name,
    description: row ? row.Description : model.description, engine: ENGINE,
    task: isClass ? 'classification' : 'regression', target, features: {columns: columnsOf(model.inputs)},
    options: reference});
}

/** Returns the names of the models whose write was refused: one refusal does not stop the others. */
async function write({register, remove}: SyncDiff): Promise<string[]> {
  const failed: string[] = [];
  for (const registration of register)
    await registerModel(registration).catch(() => failed.push(registration.model.name));
  for (const row of remove)
    await grok.functions.call('Forge:deleteModel', {model: row.Id}).catch(() => failed.push(row.Name));
  return failed;
}

function messageOf(listed: ListedModel[], failed: string[]): string {
  const messages = listed.filter(isFailedConnection).map((m) => `${m.connection} was not refreshed: ${m.error}`);
  if (failed.length > 0) {
    messages.push(`${failed.length} MLflow ${failed.length === 1 ? 'model was' : 'models were'} not updated (${
      failed[0]}${failed.length === 1 ? '' : ', ...'}).`);
  }
  return messages.join(' ');
}

/** Returns what was not refreshed. */
export async function syncModels(): Promise<string> {
  const listed = await listedModels();
  const failed = await write(syncDiff(listed, await ownRows()));
  return messageOf(listed, failed);
}

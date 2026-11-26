import * as DG from 'datagrok-api/dg';
import {EngineRegistry} from '../engines/engine-registry';
import {ForgeError} from '../forge-error';
import {forgeDb, ModelInsert, ModelRow} from '../generated/db';
import {modelNameOf, tagsOf, tagsText} from '../storage/model-fields';
import {findModelRow, isExternal, modelsChanged} from '../storage/model-store';
import {engineNamed} from '../training/train-api';
import {FeaturesSchema, TargetSchema} from '../training/train-model';
import {modelTaskOf} from './list-models';

export interface ModelInfo { name: string; description: string; tags: string[] }

/** The inputs of `Forge:updateModel`; a missing optional input is null or undefined. */
export interface UpdateModelInput {
  model: string;
  name?: string | null;
  description?: string | null;
  tags?: string | null;
}

export interface UpdatedModel { id: string; name: string; description: string; tags: string[]; version: number }

/** The body of `Forge:updateModel`: writes the given name (trimmed, not blank), description and tags (their text
 * replaces the current ones) of the visible model with the id or the unique name [UpdateModelInput.model], without a
 * version check (the last write wins). Returns the model's fields after the write. */
export async function runUpdateModel(input: UpdateModelInput): Promise<UpdatedModel> {
  const info: Partial<ModelInfo> = {};
  if (typeof input.name === 'string')
    info.name = modelNameOf(input.name);
  if (typeof input.description === 'string')
    info.description = input.description;
  if (typeof input.tags === 'string')
    info.tags = tagsOf(input.tags);
  if (Object.keys(info).length === 0)
    throw new ForgeError('Give a name, a description or tags to change.');
  const row = await findModelRow(input.model, ['name', 'description', 'tags']);
  const version = await updateModelInfo(row.id, undefined, info);
  return {id: row.id, name: info.name ?? row.name, description: info.description ?? row.description ?? '',
    tags: info.tags ?? tagsOf(row.tags), version};
}

/** Writes the given fields of the model, the tags as their column text (`tagsText`); with [version], a row changed
 * since is refused with the platform's `DomainVersionConflictError`. Returns the row's new version. */
export async function updateModelInfo(id: string, version: number | undefined,
  info: Partial<ModelInfo>): Promise<number> {
  const {tags, ...fields} = info;
  const values: Partial<ModelRow> = tags === undefined ? fields : {...fields, tags: tagsText(tags)};
  const result = await forgeDb.models.update(id, values, version === undefined ? undefined : {version});
  modelsChanged.next();
  return result.version;
}

/** The inputs of `Forge:registerModel`; a missing optional input is null or undefined. */
export interface RegisterModelInput {
  id?: string | null;
  name: string;
  description?: string | null;
  engine: string;
  task: string;
  target: TargetSchema;
  features: FeaturesSchema;
  options: {[key: string]: unknown};
}

/** The body of `Forge:registerModel`: writes the caller's own row of a model without a file and returns its id. */
export async function runRegisterModel(input: RegisterModelInput): Promise<string> {
  const engine = engineNamed(EngineRegistry.discover(), input.engine);
  const fields: ModelInsert = {name: input.name, description: input.description ?? '', engine_name: engine.name,
    engine_namespace: engine.namespace, engine_kind: 'external', task: modelTaskOf(input.task),
    target_name: input.target.name, target: input.target, features: input.features,
    feature_count: input.features.columns.length, options: input.options, storage_mode: 'none'};
  let id = input.id;
  if (id) {
    // A row the caller cannot see is not found.
    const row: ModelRow | null = await forgeDb.models.get(id);
    if (row?.author_id !== DG.User.current().id || !isExternal(row))
      throw new ForgeError(`The model '${id}' is not a model you registered.`);
    await forgeDb.models.update(id, fields);
  } else
    id = (await forgeDb.models.insert(fields))[0].id;
  modelsChanged.next();
  return id;
}

import {forgeDb, ModelRow} from '../generated/db';
import {tagsText} from '../storage/model-fields';
import {modelsChanged} from '../storage/model-store';

export interface ModelInfo { name: string; description: string; tags: string[] }

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

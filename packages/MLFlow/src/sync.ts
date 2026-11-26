/** A row of Core's `MlflowModels`: a published model, or with [error] a model that cannot be used or, without a
 * [name], a connection that could not be listed. */
export interface ListedModel {
  connectionId: string;
  connection: string;
  name: string;
  version: string;
  description: string;
  /** JSON arrays of `{name, type, semType?}`. */
  inputs: string;
  outputs: string;
  error: string;
}

/** The caller's own Forge row of a model, as `Forge:listModels` gives it. */
export interface OwnRow {
  id: string;
  name: string;
  description: string;
  target: string;
  /** The feature names, joined by ", ". */
  features: string;
  options: {connectionId?: string; name?: string; version?: string};
}

export interface Column { name: string; type: string; semType?: string }

export function columnsOf(json: string): Column[] {
  try {
    const columns: unknown = JSON.parse(json);
    return Array.isArray(columns) ? columns : [];
  } catch {
    return [];
  }
}

/** The predicted column: the model's first output. */
export function targetOf(model: ListedModel): Column {
  return columnsOf(model.outputs)[0] ?? {name: 'prediction', type: 'double'};
}

/** What a sync writes: the models to register (with [row], the row to update: its name and description are the
 * user's and stay) and the rows to delete. A row is deleted when its model is no longer listed, unless its connection
 * could not be listed; a second row of a model always. A model listed with an error is not registered, and its row
 * stays. */
export function syncDiff(listed: ListedModel[], rows: OwnRow[]): {register: {model: ListedModel; row?: OwnRow}[];
  remove: OwnRow[]} {
  const keyOf = (connectionId?: string, name?: string) => `${connectionId}/${name}`;
  const failed = new Set(listed.filter((m) => m.error !== '' && m.name === '').map((m) => m.connectionId));
  const byKey = new Map<string, OwnRow>();
  const remove: OwnRow[] = [];
  for (const row of rows) {
    const key = keyOf(row.options.connectionId, row.options.name);
    if (byKey.has(key))
      remove.push(row);
    else
      byKey.set(key, row);
  }
  const register: {model: ListedModel; row?: OwnRow}[] = [];
  for (const model of listed.filter((m) => m.name !== '')) {
    const key = keyOf(model.connectionId, model.name);
    const row = byKey.get(key);
    byKey.delete(key);
    if (model.error !== '')
      continue;
    // A registered version never changes, so the feature names are enough to see a new signature.
    const isSame = row !== undefined && row.options.version === model.version &&
      row.target === targetOf(model).name && row.features === columnsOf(model.inputs).map((c) => c.name).join(', ');
    if (!isSame)
      register.push({model, row});
  }
  for (const row of byKey.values()) {
    if (!failed.has(row.options.connectionId ?? ''))
      remove.push(row);
  }
  return {register, remove};
}

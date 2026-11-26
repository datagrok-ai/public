import * as DG from 'datagrok-api/dg';

export const ENGINE = 'MLflow';

/** MLflow's type names as Datagrok's; any other type is a string. */
export const TYPES: {[type: string]: DG.ColumnType} = {boolean: DG.COLUMN_TYPE.BOOL, integer: DG.COLUMN_TYPE.INT,
  long: DG.COLUMN_TYPE.INT, int: DG.COLUMN_TYPE.INT, float: DG.COLUMN_TYPE.FLOAT, double: DG.COLUMN_TYPE.FLOAT,
  datetime: DG.COLUMN_TYPE.DATE_TIME};

export interface Column { name: string; type: DG.ColumnType; semType?: string }

/** A row of Core's `mlflowModels` (an empty cell is absent): a model, or without [name] a connection that could not
 * be listed. */
export interface ListedModel {
  connection_id: string;
  connection: string;
  name?: string;
  version?: string;
  description?: string;
  /** JSON arrays of {@link Column}. */
  inputs?: string;
  outputs?: string;
  error?: string;
}

/** What a Forge row keeps in `options`; [type] is the declared type of the model's first output. */
export interface ModelReference { connectionId: string; name: string; version: string; type: DG.ColumnType }

/** A row of `Forge:listModels` (an empty cell is absent). */
export interface OwnRow {
  Id: string;
  Name: string;
  Description?: string;
  Target?: string;
  /** The feature names, joined by ", ". */
  Features?: string;
  /** A {@link ModelReference} as json. */
  Options: string;
}

export type NamedModel = ListedModel & {name: string};

export interface Registration {
  model: NamedModel;
  /** The row to update; without it a new one is written. */
  row?: OwnRow;
}

export interface SyncDiff { register: Registration[]; remove: OwnRow[] }

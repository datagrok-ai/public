import * as DG from 'datagrok-api/dg';
import {featureSchemasOf, targetSchemaOf} from '../apply/apply-model';
import {ApplicationStatus, ModelStorageMode, ModelTask} from '../generated/db';
import {metricFindings} from '../metrics/metrics';
import {PreparationOptions, preparationOptionsOf} from '../preparation/preparation-options';
import {datasetRefCaption, storedDatasetRef} from '../storage/dataset-ref';
import {tagsOf} from '../storage/model-fields';
import {findModelRow, ownBlob} from '../storage/model-store';
import {MetricsRecord, metricsRecordOf} from '../training/train-model';
import {applicableTables} from './applicable-tables';
import {loginsOf} from './list-models';
import {modelActivity} from './model-activity';

/** Every `model` column but `legacy_id` (the system columns come with any query). */
const GET_MODEL_COLUMNS = ['name', 'description', 'tags', 'engine_name', 'engine_namespace', 'engine_kind',
  'task', 'target_name', 'target', 'features', 'feature_count', 'options', 'hyperparameters', 'metrics', 'seed',
  'splitting', 'storage_mode', 'dataset_name', 'row_count', 'dataset_ref', 'dataset_table_id', 'dataset_fingerprint',
  'blob', 'has_training_rows'] as const;

/** What `Forge:getModel` returns: the Details, Performance and Activity panes of a model as one object. */
export interface ModelDetails {
  id: string;
  version: number;
  name: string;
  description: string;
  tags: string[];
  method: string;
  package: string;
  task: ModelTask;
  target: string;
  targetType: string | null;
  /** A classifier's categories; null without them. */
  classes: string[] | null;
  /** In training order. */
  features: string[];
  featureTypes: {[name: string]: string};
  hyperparameters: {[name: string]: unknown};
  preparation: PreparationOptions;
  metrics: MetricsRecord | null;
  findings: {text: string; detail: string}[];
  seed: number | null;
  splitting: {[key: string]: unknown} | null;
  storage: {mode: ModelStorageMode; source?: string; tableId?: string};
  dataset: {name: string; rows: number | null; summary: {[key: string]: unknown} | null};
  hasModelFile: boolean;
  keepsTrainingRows: boolean;
  /** The author's login; empty when it cannot be read. */
  author: string;
  created: string;
  updated: string;
  lastRun: {when: string; status: ApplicationStatus} | null;
  applications: number;
  /** The names of the [tables] the model fits. */
  applicableTo: string[];
}

/** The body of `Forge:getModel`: the visible model with the id or the unique name [idOrName], read with its
 * applications and author; [tables] are the open tables `applicableTo` checks. */
export async function runGetModel(idOrName: string, tables: DG.DataFrame[]): Promise<ModelDetails> {
  const row = await findModelRow(idOrName, GET_MODEL_COLUMNS);
  const [activity, logins] = await Promise.all([modelActivity(row.id), loginsOf([row.author_id])]);
  const schemas = featureSchemasOf(row.features) ?? [];
  const target = targetSchemaOf(row.target);
  const metrics = metricsRecordOf(row.metrics);
  const storage: ModelDetails['storage'] = {mode: row.storage_mode};
  const ref = storedDatasetRef(row.dataset_ref);
  if (ref !== null)
    storage.source = datasetRefCaption(ref);
  if (row.dataset_table_id)
    storage.tableId = row.dataset_table_id;
  return {
    id: row.id,
    version: row.version,
    name: row.name,
    description: row.description ?? '',
    tags: tagsOf(row.tags),
    method: row.engine_name,
    package: row.engine_namespace,
    task: row.task,
    target: row.target_name,
    targetType: target?.type ?? null,
    classes: target?.categories ?? null,
    features: schemas.map((f) => f.name),
    featureTypes: Object.fromEntries(schemas.map((f) => [f.name, f.type])),
    hyperparameters: row.hyperparameters ?? {},
    preparation: preparationOptionsOf(row.options),
    metrics,
    findings: metrics === null ? [] : metricFindings(metrics).map(({text, detail}) => ({text, detail})),
    seed: row.seed ?? null,
    splitting: row.splitting ?? null,
    storage,
    dataset: {name: row.dataset_name ?? '', rows: row.row_count ?? null, summary: row.dataset_fingerprint ?? null},
    hasModelFile: ownBlob(row.blob) !== null,
    keepsTrainingRows: row.has_training_rows === true,
    author: logins.get(row.author_id) ?? '',
    created: row.created_on.toISOString(),
    updated: row.updated_on.toISOString(),
    lastRun: activity.lastRun === undefined ? null :
      {when: activity.lastRun.when.toISOString(), status: activity.lastRun.status},
    applications: activity.count,
    applicableTo: applicableTables(row, tables).map((t) => t.name),
  };
}

import {Engine, hasTrainingRows} from '../engines/engine';
import {ModelInsert, TrainingRunInsert, TrainingRunStatus} from '../generated/db';
import {PreparationOptions} from '../preparation/preparation-options';
import {FeaturesSchema, MetricsRecord, Splitting, TargetSchema, TrainingRequest, TrainingResult, trainingSetupOf}
  from '../training/train-model';
import {DatasetFingerprint} from './dataset-fingerprint';
import {DatasetRef} from './dataset-ref';

type JsonColumn = 'target' | 'features' | 'options' | 'metrics' | 'splitting' | 'dataset_fingerprint' | 'dataset_ref';

export type ModelFields = Omit<ModelInsert, JsonColumn> & {
  target: TargetSchema;
  features: FeaturesSchema;
  options: PreparationOptions;
  metrics: MetricsRecord;
  splitting: Splitting;
  dataset_fingerprint: DatasetFingerprint;
  dataset_ref?: DatasetRef;
};

/** What a model keeps of its training data: only the fingerprint, a link to the source, or an uploaded copy. */
export type ModelStorage = {mode: 'none'} | {mode: 'reference'; ref: DatasetRef} | {mode: 'copy'; tableId: string};

export type TrainingRunRecord = Omit<TrainingRunInsert, JsonColumn> & {
  features: FeaturesSchema;
  options: PreparationOptions;
  metrics?: MetricsRecord;
  splitting: Splitting;
  dataset_fingerprint: DatasetFingerprint;
};

export function modelFieldsOf(input: {name: string; description: string; tags: string[]; engine: Engine;
  datasetName: string; result: TrainingResult; fingerprint: DatasetFingerprint; storage: ModelStorage}): ModelFields {
  const {engine, result, storage} = input;
  return {
    name: input.name,
    description: input.description,
    tags: tagsText(input.tags),
    ...engineFieldsOf(engine),
    task: result.task,
    target_name: result.target.name,
    target: result.target,
    features: result.features,
    feature_count: result.features.columns.length,
    options: result.options,
    hyperparameters: result.hyperparameters,
    metrics: result.metrics,
    seed: result.seed,
    splitting: result.splitting,
    storage_mode: storage.mode,
    ...(storage.mode === 'reference' ? {dataset_ref: storage.ref} : {}),
    ...(storage.mode === 'copy' ? {dataset_table_id: storage.tableId} : {}),
    dataset_name: input.datasetName,
    row_count: result.rowCount,
    dataset_fingerprint: input.fingerprint,
    has_training_rows: hasTrainingRows(engine),
  };
}

export function trainingRunOf(input: {request: TrainingRequest; datasetName: string; fingerprint: DatasetFingerprint;
  status: TrainingRunStatus; startedOn: string; durationMs: number; metrics?: MetricsRecord;
  error?: string}): TrainingRunRecord {
  const {request} = input;
  const setup = trainingSetupOf(request);
  return {
    ...engineFieldsOf(request.engine),
    task: setup.task,
    target_name: request.target.name,
    features: setup.features,
    options: setup.options,
    hyperparameters: request.hyperparameters,
    seed: request.seed,
    splitting: setup.splitting,
    metrics: input.metrics,
    dataset_name: input.datasetName,
    row_count: request.target.length,
    dataset_fingerprint: input.fingerprint,
    status: input.status,
    error: input.error,
    started_on: input.startedOn,
    duration_ms: input.durationMs,
  };
}

/** The tags as stored: trimmed, without empty ones and repeats, in the given order. */
export function normalizedTags(tags: string[]): string[] {
  return [...new Set(tags.map((t) => t.trim()).filter((t) => t !== ''))];
}

/** The `tags` column text of [tags]: joined by ", "; empty without tags, which the server stores as null. */
export function tagsText(tags: string[]): string {
  return normalizedTags(tags).join(', ');
}

/** The tags of a `tags` column value; every comma separates two tags. */
export function tagsOf(text: string | null | undefined): string[] {
  return normalizedTags((text ?? '').split(','));
}

function engineFieldsOf(engine: Engine): Pick<ModelInsert, 'engine_name' | 'engine_namespace' | 'engine_kind'> {
  return {engine_name: engine.name, engine_namespace: engine.namespace, engine_kind: engine.kind};
}

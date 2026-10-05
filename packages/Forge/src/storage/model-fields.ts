import {Engine} from '../engines/engine';
import {ModelInsert, TrainingRunInsert, TrainingRunStatus} from '../generated/db';
import {PreparationOptions} from '../preparation/preparation-options';
import {FeaturesSchema, MetricsRecord, Splitting, TargetSchema, TrainingRequest, TrainingResult, trainingSetupOf}
  from '../training/train-model';
import {DatasetFingerprint} from './dataset-fingerprint';

type JsonColumn = 'target' | 'features' | 'options' | 'metrics' | 'splitting' | 'dataset_fingerprint';

export type ModelFields = Omit<ModelInsert, JsonColumn> & {
  target: TargetSchema;
  features: FeaturesSchema;
  options: PreparationOptions;
  metrics: MetricsRecord;
  splitting: Splitting;
  dataset_fingerprint: DatasetFingerprint;
};

export type TrainingRunRecord = Omit<TrainingRunInsert, JsonColumn> & {
  features: FeaturesSchema;
  options: PreparationOptions;
  metrics?: MetricsRecord;
  splitting: Splitting;
  dataset_fingerprint: DatasetFingerprint;
};

export function modelFieldsOf(input: {name: string; description: string; engine: Engine; datasetName: string;
  result: TrainingResult; fingerprint: DatasetFingerprint}): ModelFields {
  const {engine, result} = input;
  return {
    name: input.name,
    description: input.description,
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
    storage_mode: 'none',
    dataset_name: input.datasetName,
    row_count: result.rowCount,
    dataset_fingerprint: input.fingerprint,
    has_training_rows: false,
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

function engineFieldsOf(engine: Engine): Pick<ModelInsert, 'engine_name' | 'engine_namespace' | 'engine_kind'> {
  return {engine_name: engine.name, engine_namespace: engine.namespace, engine_kind: engine.kind};
}

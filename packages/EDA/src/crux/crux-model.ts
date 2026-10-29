import * as DG from 'datagrok-api/dg';

import {cruxMl} from './crux-loader';
import {floatColumn, listColumns, stringValues, toArrowIpc} from './crux-bridge';
import {LEGACY_MODEL_ERROR, PREDICTION, blobId, featureColumns, runCrux} from './model-common';
import {legacyPredict} from '../legacy/legacy-models';

export {LEGACY_MODEL_ERROR, runCrux};

const PROBA_PREFIX = 'proba_';

const FORMAT = 'crux-model-v1';
const MANIFEST_KEY = 'crux.model';
const FEATURES_KEY = 'eda.features';

export type CruxTask = 'regression' | 'classification';

export interface CruxFitParams {
  estimator?: Record<string, number | string | boolean | number[]>;
  preprocess?: 'standard' | 'minmax' | 'maxabs' | 'robust';
  seed?: number;
}

export interface CruxManifest {
  format: string;
  task: CruxTask;
  model: string;
  classes?: (string | boolean | number)[];
  n_features: number;
  pipeline?: {step: string}[];
  params: Record<string, unknown>;
  training: Record<string, unknown>;
  eda?: {features: string[]};
}

function safetensorsHeader(bytes: Uint8Array): {header: any, data: Uint8Array} | null {
  if (bytes.length < 8)
    return null;
  const n = Number(new DataView(bytes.buffer, bytes.byteOffset, 8).getBigUint64(0, true));
  if (n <= 0 || n > bytes.length - 8)
    return null;
  try {
    return {header: JSON.parse(new TextDecoder().decode(bytes.subarray(8, 8 + n))), data: bytes.subarray(8 + n)};
  } catch {
    return null;
  }
}

export function isCruxModel(bytes: Uint8Array): boolean {
  const md = safetensorsHeader(bytes)?.header?.__metadata__;
  return md?.format === FORMAT && typeof md[MANIFEST_KEY] === 'string';
}

/** The manifest of a crux-model-v1 blob; any other blob throws {@link LEGACY_MODEL_ERROR} (`cruxPredict` routes
 * such blobs to `legacyPredict`) */
export function readCruxManifest(bytes: Uint8Array): CruxManifest {
  if (!isCruxModel(bytes))
    throw new Error(LEGACY_MODEL_ERROR);
  const md = safetensorsHeader(bytes)!.header.__metadata__;
  const manifest = JSON.parse(md[MANIFEST_KEY]) as CruxManifest;
  if (md[FEATURES_KEY] !== undefined)
    manifest.eda = {features: JSON.parse(md[FEATURES_KEY])};
  return manifest;
}

/** Adds the DG feature names to the safetensors `__metadata__`; the blob stays a crux-model-v1 file */
function withEdaMetadata(bytes: Uint8Array, features: string[]): Uint8Array {
  const {header, data} = safetensorsHeader(bytes)!;
  header.__metadata__[FEATURES_KEY] = JSON.stringify(features);
  const json = new TextEncoder().encode(JSON.stringify(header));
  const len = Math.ceil(json.length / 8) * 8;
  const out = new Uint8Array(8 + len + data.length);
  new DataView(out.buffer).setBigUint64(0, BigInt(len), true);
  out.fill(0x20, 8, 8 + len);
  out.set(json, 8);
  out.set(data, 8 + len);
  return out;
}

/** Trains a Crux model on all columns of `df`; returns a crux-model-v1 blob */
export async function cruxFit(task: CruxTask, model: string, df: DG.DataFrame, predictColumn: DG.Column,
  params: CruxFitParams = {}): Promise<Uint8Array> {
  const features = df.columns.toList();
  const seed = params.seed === undefined ? {} : {random_state: params.seed};
  const callParams = {task, model, featurizer: 'vector', cv: 0, preprocess: params.preprocess,
    estimator: JSON.stringify({...params.estimator, ...seed}), ...seed};

  const out = await runCrux('ml.fit', toArrowIpc(features, predictColumn), callParams);
  return withEdaMetadata(out.getChild('model')!.get(0) as Uint8Array, features.map((c) => c.name));
}

/** The columns of `df` the model was trained on, matched by the names saved at training */
export function modelFeatures(df: DG.DataFrame, manifest: CruxManifest): DG.Column[] {
  const names = manifest.eda?.features;
  if (names === undefined)
    throw new Error('Crux: the model has no feature names (not trained by EDA)');
  return featureColumns(df, names);
}

/** A classifier whose `ml.predict` fills `probas` (not `svc`), with column names that differ beyond case (a DataFrame
 * refuses two names that differ only in case) */
function hasProbabilities(manifest: CruxManifest): boolean {
  const classes = (manifest.classes ?? []).map((c) => String(c));
  return manifest.model !== 'svc' && new Set(classes.map((c) => c.toLowerCase())).size === classes.length;
}

/** Applies a crux-model-v1 blob or a pre-Crux EDA model; features are matched by the names saved at training */
export async function cruxPredict(df: DG.DataFrame, bytes: Uint8Array): Promise<DG.DataFrame> {
  if (!isCruxModel(bytes))
    return legacyPredict(df, bytes);
  const manifest = readCruxManifest(bytes);
  const columns = modelFeatures(df, manifest);

  const out = await runCrux('ml.predict', toArrowIpc(columns), {}, {name: 'model', id: blobId(bytes), bytes});

  if (manifest.task === 'regression')
    return DG.DataFrame.fromColumns([floatColumn(out, 'value', PREDICTION)]);

  const classes = manifest.classes ?? [];
  const probabilities = hasProbabilities(manifest) ?
    listColumns(out, 'probas', classes.map((c) => `${PROBA_PREFIX}${c}`)) : [];
  const labels = stringValues(out, 'class');
  const all = (type: string) => classes.length > 0 && classes.every((c) => typeof c === type);
  const result = (type: DG.ColumnType, values: any[]) =>
    DG.DataFrame.fromColumns([DG.Column.fromList(type, PREDICTION, values), ...probabilities]);
  if (all('boolean'))
    return result(DG.COLUMN_TYPE.BOOL, labels.map((v) => v === 'true'));
  if (all('number')) {
    const type = classes.every((c) => Number.isInteger(c)) ? DG.COLUMN_TYPE.INT : DG.COLUMN_TYPE.FLOAT;
    return result(type, labels.map((v) => v === null ? null : Number(v)));
  }
  return result(DG.COLUMN_TYPE.STRING, labels);
}

export function cruxTerminate(): void {
  cruxMl.terminate();
}

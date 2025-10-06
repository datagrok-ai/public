import * as DG from 'datagrok-api/dg';

/* eslint-disable no-unused-vars */
export enum ErrorHandling {
  REJECT_ALL = 'reject_all',
  REJECT_ROW = 'reject_row'
}

export enum Scope {
  COMPOUNDS = 'compounds',
  BATCHES = 'batches',
  ASSAYS = 'assays',
  ASSAY_RUNS = 'assay_runs',
  ASSAY_RESULTS = 'assay_results'
}

export const ErrorHandlingLabels: Record<string, ErrorHandling> = {
  'Reject all': ErrorHandling.REJECT_ALL,
  'Reject invalid': ErrorHandling.REJECT_ROW,
};

export const ScopeLabels: Record<string, Scope> = {
  'Compounds': Scope.COMPOUNDS,
  'Batches': Scope.BATCHES,
  'Assays': Scope.ASSAYS,
  'Assay Runs': Scope.ASSAY_RUNS,
  'Assay Results': Scope.ASSAY_RESULTS,
};

export const ScopeLabelsReduced: Record<string, Scope> = {
  'Compounds': Scope.COMPOUNDS,
  'Batches': Scope.BATCHES,
};

export enum ResultOutput {
  CSV = 'csv',
  JSON = 'json'
}

export const scopeToUrl: { [key: string]: string } = {
  compounds: '/v1/compounds/',
  batches: '/v1/batches/',
  assays: '/v1/assays',
  assay_runs: '/v1/assay_runs/',
  assay_results: '/v1/assay_results/',
};

export type MolTrackProp = {
  name: string;
  value_type: string;
  entity_type?: string;
  description?: string;
  pattern?: string;
  friendly_name?: string;
};

export const GITHUB_BASE_URL =
  'https://raw.githubusercontent.com/datagrok-ai/mol-track/main/data/black/';


export const EXCLUDE_SEARCH_FIELDS = ['id', 'molregno', 'batch_regno', 'original_molfile', 'hash_mol',
  'hash_tautomer', 'created_by', 'updated_by', 'deleted_by', 'hash_canonical_smiles', 'hash_no_stereo_smiles',
  'hash_no_stereo_tautomer'];

export const EXCLUDE_SEARCH_OUTPUT_FIELDS = ['id', 'molregno', 'batch_regno', 'hash_mol', 'hash_tautomer', 'uuid',
  'hash_canonical_smiles', 'hash_no_stereo_smiles', 'hash_no_stereo_tautomer'];

export const MOL_COL_NAME = 'canonical_smiles';
export const CORPORATE_COMPOUND_ID_COL_NAME = 'corporate_compound_id';
export const STRUCTURE_FIELDS = [MOL_COL_NAME, 'original_molfile'];
export const STRUCTURE_SEARCH_FIELD = 'structure';

export const STRING_AGGREGATIONS = ['CONCAT ALL', 'CONCAT UNIQUE', 'LONGEST', 'SHORTEST', 'MOST FREQUENT'];
export const NUMERIC_AGGREGATIONS = ['FIRST', 'COUNT', 'VALUES', 'UNIQUE', 'NULLS', 'MIN', 'MAX', 'SUM',
  'MED', 'AVG', 'STDEV', 'VARIANCE', 'Q1', 'Q2', 'Q3'];
export const PROP_NUM_TYPES = [DG.TYPE.BIG_INT, DG.TYPE.FLOAT, DG.TYPE.INT, DG.TYPE.NUM, DG.TYPE.QNUM];

export const MOLTRACK_ENDPOINT = 'moltrackEndpoint';
export const MOLTRACK_ENTITY_TYPE = 'moltrackEntityType';
export const MOLTRACK_ENTITY_LEVEL = 'moltrackEntityLevel';
export const MOLTRACK_IS_STATIC_FIELD = 'moltrackIsStaticField';
export const PROPERTIES = 'properties';

export const SEARCH_NODE = 'Search';
export const SAVED_SEARCHES_NODE = 'Saved Searches';

export const MOLTRACK_APP_PATH: string = 'apps/MolTrack';

export const MOLTRACK_REQUEST_TITLE_UPDATE = 'moltrack-updateTitleRequested';
export const MOLTRACK_MAPPING_VALIDATION_CHANGED = 'moltrack-mappingValidationChanged';

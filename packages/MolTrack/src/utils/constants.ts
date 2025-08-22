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
  description?: string;
  pattern?: string;
};

export const GITHUB_BASE_URL =
  'https://raw.githubusercontent.com/datagrok-ai/mol-track/main/data/black/';

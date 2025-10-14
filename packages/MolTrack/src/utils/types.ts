import { ComplexCondition, Operators } from '@datagrok-libraries/utils/src/query-builder/query-builder';
import { Scope } from './constants';
import * as DG from 'datagrok-api/dg';

/* eslint-disable no-unused-vars */

// Enums for property classification and entity types
export enum MolTrackPropertyClass {
  DECLARED = 'DECLARED',
  CALCULATED = 'CALCULATED',
  MEASURED = 'MEASURED',
  PREDICTED = 'PREDICTED'
}

export enum MolTrackEntityType {
  COMPOUND = 'COMPOUND',
  BATCH = 'BATCH',
  ASSAY = 'ASSAY',
  ASSAY_RUN = 'ASSAY_RUN',
  ASSAY_RESULT = 'ASSAY_RESULT',
}

// MolTrack Property interface
export interface MolTrackProperty {
  description: string;
  entity_type: MolTrackEntityType;
  name: string;
  friendly_name: string;
  pattern?: string;
  property_class: MolTrackPropertyClass;
  semantic_type_id: number;
  value_type: 'string' | 'int' | 'double' | 'datetime' | 'uuid' | 'boolean';
  unit?: string;
  semantic_type?: MolTrackSemanticType;
  min?: number;
  max?: number;
  choices?: any;
  validators?: any;
  value_datetime?: any;
  value_uuid?: string;
  value_num?: number;
  value_string?: string;
}

export interface MolTrackSemanticType {
  name?: string;
  description?: string;
}

// Simple condition structure
export interface MolTrackSimpleCondition {
  field: string;
  operator: string;
  value: any;
  threshold?: number | null; // Only used with IS SIMILAR operator
}

// Complex condition structure
export interface MolTrackComplexCondition {
  operator: Operators.Logical;
  conditions: MolTrackFilter[];
}

// Union type for any condition
export type MolTrackFilter = MolTrackSimpleCondition | MolTrackComplexCondition;

export interface MolTrackSearchAggregation {
   field: string;
   operation: string;
}
export interface MolTrackSearchQuery {
  level: string;
  output: string[];
  filter?: MolTrackFilter;
  aggregations?: MolTrackSearchAggregation[];
}

export interface MolTrackSearch {
  outputCols: {name: string, type: string}[];
  condition: ComplexCondition;
  aggregations?: MolTrackSearchAggregation[];
}

export interface MolTrackSearchResponse {
    status: string,
    data: Record<string, any>[],
    total_count: number,
    level: string,
    columns: string[],
}

export type MolTrackSearchHistoryItem = {
    date: DG.TYPE.DATE_TIME,
    value: string,
}

export const molTrackSearchMapping = {
  [MolTrackEntityType.COMPOUND]: {level: Scope.COMPOUNDS, searchEndpoint: 'compounds'},
  [MolTrackEntityType.BATCH]: {level: Scope.BATCHES, searchEndpoint: 'batches'},
  [MolTrackEntityType.ASSAY]: {level: Scope.ASSAYS, searchEndpoint: 'assays'},
  [MolTrackEntityType.ASSAY_RUN]: {level: Scope.ASSAY_RUNS, searchEndpoint: 'assay-runs'},
  [MolTrackEntityType.ASSAY_RESULT]: {level: Scope.ASSAY_RESULTS, searchEndpoint: 'assay-results'},
};



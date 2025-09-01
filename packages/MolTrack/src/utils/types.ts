import { Operators } from '@datagrok-libraries/utils/src/query-builder/query-builder';
import { Scope } from './constants';


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
  semantic_type?: MolTrackSemanticType
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

// Complex condition structure
export interface MolTrackSearchQuery {
  level: string;
  output: string[];
  filter: MolTrackFilter;
}

export interface MolTrackSearchResponse {
    status: string,
    data: Record<string, any>[],
    total_count: number,
    level: string,
    columns: string[],
}

export interface MolTrackSearchMapping {
    level: Scope,
    searchEndpoint: string,
    propEntityType: MolTrackEntityType
}

export const searchTypeMapping: MolTrackSearchMapping[] = [
  { level: Scope.COMPOUNDS, searchEndpoint: 'compounds', propEntityType: MolTrackEntityType.COMPOUND },
  { level: Scope.BATCHES, searchEndpoint: 'batches', propEntityType: MolTrackEntityType.BATCH },
  { level: Scope.ASSAYS, searchEndpoint: 'assays', propEntityType: MolTrackEntityType.ASSAY },
  { level: Scope.ASSAY_RUNS, searchEndpoint: 'assay-runs', propEntityType: MolTrackEntityType.ASSAY_RUN },
  { level: Scope.ASSAY_RESULTS, searchEndpoint: 'assay-results', propEntityType: MolTrackEntityType.ASSAY_RESULT },
];


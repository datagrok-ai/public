import {KnimeJobState} from './constants';

export type KnimeParamType = 'table' | 'string' | 'int' | 'double' | 'boolean' | 'file' | 'json';

export interface KnimeInputParam {
  name: string;
  type: KnimeParamType;
  required: boolean;
  description?: string;
  /** Default value from the OpenAPI spec. */
  defaultValue?: any;
  /** Parent parameter name for grouped inputs (variable/json sub-properties). */
  group?: string;
  /** Description of the parent group (from the OpenAPI spec). */
  groupDescription?: string;
  /** Table column schema from the OpenAPI spec example (for table inputs). */
  tableSpec?: {name: string; type: string}[];
}

export interface KnimeWorkflowSpec {
  inputs: KnimeInputParam[];
  hasFileOutputs: boolean;
}

export interface KnimeTableInput {
  'table-spec': {name: string; type: string}[];
  'table-data': any[][];
}

export interface KnimeExecutionInput {
  [paramName: string]: any;
}

export interface KnimeExecutionResult {
  outputs: {[paramName: string]: any};
  outputTables?: any[];
  outputResources?: {[resourceId: string]: string};
  /** Pre-fetched output resources (fetched immediately to avoid job discard race). */
  fetchedResources?: {name: string; resource: KnimeOutputResource}[];
  /** Job ID from the execution response (available even for sync REST execution). */
  jobId?: string;
  dataAppUrl?: string;
  error?: string;
}

export interface KnimeJobStatus {
  id: string;
  state: KnimeJobState;
  message?: string;
  progress?: number;
  /** Raw API response data (used internally to avoid re-fetching). */
  _rawData?: any;
}

export interface KnimeOutputResource {
  contentType: string;
  /** Parsed JSON data (when content-type is application/json). */
  json?: any;
  /** Text content (when content-type is text/*). */
  text?: string;
  /** Raw binary data (when content is not JSON or text). */
  blob?: Blob;
}

export interface KnimeDeployment {
  id: string;
  name: string;
  type: string;
  state?: string;
  workflowPath?: string;
}

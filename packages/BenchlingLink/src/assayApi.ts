import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { dataFrameFromObjects } from './utils';

// Types for Assay Results and Assay Runs API
export interface AssayResult {
  archiveRecord?: any;
  createdAt?: string;
  creator?: any;
  entryId?: string | null;
  fieldValidation?: any;
  fields?: any;
  id: string;
  isReviewed?: boolean;
  modifiedAt?: string;
  projectId?: string | null;
  schema?: any;
  validationComment?: string;
  validationStatus?: string;
}

export interface AssayRun {
  apiURL?: string;
  archiveRecord?: any;
  createdAt?: string;
  creator?: any;
  entryId?: string | null;
  equipmentId?: string | null;
  fields?: any;
  id: string;
  isReviewed?: boolean;
  projectId?: string | null;
  schema?: any;
  validationComment?: string;
  validationStatus?: string;
}

export interface AssayResultsQueryParams {
  schemaId?: string;
  createdAt_lt?: string;
  createdAt_gt?: string;
  createdAt_lte?: string;
  createdAt_gte?: string;
  minCreatedTime?: number;
  maxCreatedTime?: number;
  sort?: string;
  nextToken?: string;
  pageSize?: number;
  entityIds?: string;
  storageIds?: string;
  assayRunIds?: string;
  automationOutputProcessorId?: string;
  ids?: string;
  modifiedAt_lt?: string;
  modifiedAt_gt?: string;
  modifiedAt_lte?: string;
  modifiedAt_gte?: string;
  archiveReason?: string;
}

export interface AssayRunsQueryParams {
  schemaId: string;
  minCreatedTime?: number;
  maxCreatedTime?: number;
  nextToken?: string;
  pageSize?: number;
  ids?: string;
}

export async function queryAssayResults(params: AssayResultsQueryParams = {}): Promise<DG.DataFrame> {
  // const token = 'YOUR_BENCHLING_API_TOKEN';
  // const query = buildQueryString(params);
  // const url = `https://benchling.com/api/v2/assay-results${query ? '?' + query : ''}`;
  // const response = await grok.dapi.fetchProxy(url, {
  //   headers: {
  //     'Authorization': `Bearer ${token}`,
  //     'Accept': 'application/json',
  //   },
  // });
  // if (!response.ok)
  //   throw new Error(`Benchling API error: ${response.statusText}`);
  // const data = await response.json();
  // const df = dataFrameFromObjects(data.assayResults ?? []) ?? DG.DataFrame.create();
  // return df;
  const df = dataFrameFromObjects(mockAssayResults.assayResults);
  return df;
}


export async function queryAssayRuns(params: AssayRunsQueryParams): Promise<DG.DataFrame> {
  // const token = 'YOUR_BENCHLING_API_TOKEN';
  // const query = buildQueryString(params);
  // const url = `https://benchling.com/api/v2/assay-runs${query ? '?' + query : ''}`;
  // const response = await grok.dapi.fetchProxy(url, {
  //   headers: {
  //     'Authorization': `Bearer ${token}`,
  //     'Accept': 'application/json',
  //   },
  // });
  // if (!response.ok)
  //   throw new Error(`Benchling API error: ${response.statusText}`);
  // const data = await response.json();
  // const df = dataFrameFromObjects(data.assayRuns ?? []) ?? DG.DataFrame.create();
  // return df;
  const df = dataFrameFromObjects(mockAssayRuns.assayRuns);
  return df;
}


export interface AssayResultCreate {
  schemaId: string;
  fields: any;
  id?: string;
  projectId?: string | null;
  fieldValidation?: any;
}

export interface AssayResultsBulkCreateRequest {
  assayResults: AssayResultCreate[];
}

export interface AssayResultsCreateResponse {
  assayResults: string[];
  errors?: any;
}

export async function postAssayResults(body: AssayResultsBulkCreateRequest): Promise<AssayResultsCreateResponse> {
  // const token = 'YOUR_BENCHLING_API_TOKEN';
  // const response = await grok.dapi.fetchProxy('https://benchling.com/api/v2/assay-results', {
  //   method: 'POST',
  //   headers: {
  //     'Authorization': `Bearer ${token}`,
  //     'Accept': 'application/json',
  //     'Content-Type': 'application/json',
  //   },
  //   body: JSON.stringify(body),
  // });
  // if (!response.ok)
  //   throw new Error(`Benchling API error: ${response.statusText}`);
  // const data = await response.json();
  // return data;
  const count = body.assayResults?.length ?? 0;
  return { assayResults: mockAssayResults.assayResults.slice(0, count).map((r) => r.id) };
}

export interface AssayRunCreate {
  schemaId: string;
  fields: any;
  id?: string;
  projectId?: string;
  validationComment?: string;
  validationStatus?: string;
}

export interface AssayRunsBulkCreateRequest {
  assayRuns: AssayRunCreate[];
}

export interface AssayRunsBulkCreateResponse {
  assayRuns: string[];
  errors?: any;
}

export async function postAssayRuns(body: AssayRunsBulkCreateRequest): Promise<AssayRunsBulkCreateResponse> {
  // const token = 'YOUR_BENCHLING_API_TOKEN';
  // const response = await grok.dapi.fetchProxy('https://benchling.com/api/v2/assay-runs', {
  //   method: 'POST',
  //   headers: {
  //     'Authorization': `Bearer ${token}`,
  //     'Accept': 'application/json',
  //     'Content-Type': 'application/json',
  //   },
  //   body: JSON.stringify(body),
  // });
  // if (!response.ok)
  //   throw new Error(`Benchling API error: ${response.statusText}`);
  // const data = await response.json();
  // return data;
  const count = body.assayRuns?.length ?? 0;
  return { assayRuns: mockAssayRuns.assayRuns.slice(0, count).map((r) => r.id) };
}

export const mockAssayResults = {
  assayResults: [
    ...Array.from({length: 100}, (_, i) => ({
      id: `ar_${String(i+1).padStart(3, '0')}`,
      createdAt: `2023-${String((i%12)+1).padStart(2, '0')}-01T12:00:00Z`,
      modifiedAt: `2023-${String((i%12)+1).padStart(2, '0')}-02T12:00:00Z`,
      isReviewed: i % 2 === 0,
      fields: { result: 10 + i, unit: 'ng/mL' },
      schema: { id: `schema_${i+1}`, name: `ResultSchema${i+1}` },
    }))
  ],
  nextToken: undefined,
};

export const mockAssayRuns = {
  assayRuns: [
    ...Array.from({length: 100}, (_, i) => ({
      id: `run_${String(i+1).padStart(3, '0')}`,
      createdAt: `2023-${String((i%12)+1).padStart(2, '0')}-01T12:00:00Z`,
      isReviewed: i % 2 === 0,
      fields: { runValue: 100 + i, status: i % 2 === 0 ? 'complete' : 'pending' },
      schema: { id: `schema_${i+1}`, name: `RunSchema${i+1}` },
      apiURL: `https://benchling.com/api/v2/assay-runs/run_${String(i+1).padStart(3, '0')}`,
    }))
  ],
  nextToken: undefined,
};
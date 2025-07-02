// Types for AA Sequences API
export interface AaAnnotation {
  color?: string;
  end: number;
  id: string;
  name: string;
  notes?: string;
  start: number;
  type: string;
}

export interface UserSummary {
  handle: string;
  id: string;
  name: string;
}

export interface AaSequence {
  aliases?: string[];
  aminoAcids: string;
  annotations?: AaAnnotation[];
  apiURL?: string;
  archiveRecord?: any;
  createdAt?: string;
  creator?: UserSummary;
  customFields?: any;
  entityRegistryId?: string | null;
  fields?: any;
  folderId?: string | null;
  id: string;
  length: number;
  modifiedAt?: string;
  name: string;
  registrationOrigin?: any;
  registryId?: string | null;
  schema?: any;
  webURL?: string;
}

export interface AaSequencesPaginatedList {
  aaSequences: AaSequence[];
  nextToken?: string;
}

import { mockAaSequences, mockAaSequence } from './aaSequencesMock';

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { dataFrameFromObjects } from './utils';

export async function getToken(): Promise<DG.DataFrame> {
  const response = await grok.dapi.fetchProxy('https://benchling.com/api/v2/token', {
    headers: {
      'Accept': 'application/json',
    },
  });
  if (!response.ok)
    throw new Error(`Benchling API error: ${response.statusText}`);
  const data = await response.json();
  const df = DG.DataFrame.fromObjects(data.aaSequences ?? []) ?? DG.DataFrame.create();
  grok.shell.addTableView(df);
  return df;
}

export function buildQueryString(params: Record<string, any>): string {
  const esc = encodeURIComponent;
  return Object.entries(params)
    .filter(([_, v]) => v !== undefined && v !== null && v !== '')
    .map(([k, v]) => `${esc(k)}=${esc(v)}`)
    .join('&');
}

export interface AASequencesQueryParams {
  pageSize?: number;
  nextToken?: string;
  sort?: string;
  createdAt?: string;
  modifiedAt?: string;
  name?: string;
  nameIncludes?: string;
  aminoAcids?: string;
  folderId?: string;
  mentionedIn?: string;
  projectId?: string;
  registryId?: string;
  schemaId?: string;
  schemaFields?: string;
  archiveReason?: string;
  mentions?: string;
  ids?: string;
  entityRegistryIds_anyOf?: string;
  names_anyOf?: string;
  names_anyOf_caseSensitive?: string;
  creatorIds?: string;
  authorIds_anyOf?: string;
  returning?: string;
}

export async function queryAASequences(params: AASequencesQueryParams = {}): Promise<DG.DataFrame> {
  // const token = 'YOUR_BENCHLING_API_TOKEN';
  // const query = buildQueryString(params);
  // const url = `https://benchling.com/api/v2/aa-sequences${query ? '?' + query : ''}`;
  // const response = await grok.dapi.fetchProxy(url, {
  //   headers: {
  //     'Authorization': `Bearer ${token}`,
  //     'Accept': 'application/json',
  //   },
  // });
  // if (!response.ok)
  //   throw new Error(`Benchling API error: ${response.statusText}`);
  // const data = await response.json();
  // const df = dataFrameFromObjects(data.aaSequences ?? []) ?? DG.DataFrame.create();
  const df = dataFrameFromObjects(mockAaSequences.aaSequences);
  return df;
}


export interface AaSequenceCreateRequest {
  name: string;
  aminoAcids: string;
  aliases?: string[];
  annotations?: AaAnnotation[];
  authorIds?: string[];
  customFields?: any;
  fields?: any;
  folderId?: string;
  schemaId?: string;
}

/**
 * Create a new AA Sequence via POST /aa-sequences
 * @param body The request body for creating an AA sequence
 * @returns The created AA sequence
 */
export async function postAASequence(body: AaSequenceCreateRequest): Promise<AaSequence> {
  // const token = 'YOUR_BENCHLING_API_TOKEN';
  // const response = await grok.dapi.fetchProxy('https://benchling.com/api/v2/aa-sequences', {
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
  // return data as AaSequence;
  return mockAaSequences.aaSequences[0];
} 
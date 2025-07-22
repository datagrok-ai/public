import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { dataFrameFromObjects } from './utils';
import { UserSummary, ArchiveRecord } from './types';

export interface DnaAnnotation {
  color?: string;
  end: number;
  id: string;
  name: string;
  notes?: string;
  start: number;
  type: string;
}

export interface DnaOligo {
  aliases?: string[];
  annotations?: DnaAnnotation[];
  apiURL?: string;
  archiveRecord?: ArchiveRecord | null;
  bases: string;
  createdAt?: string;
  creator?: UserSummary;
  customFields?: any;
  customNotation?: string | null;
  customNotationName?: string | null;
  entityRegistryId?: string | null;
  fields?: any;
  folderId?: string | null;
  helm?: string;
  id: string;
  length: number;
  modifiedAt?: string;
  name: string;
  nucleotideType: 'DNA';
  registrationOrigin?: any;
  registryId?: string | null;
  schema?: any;
  webURL?: string;
}

export interface DnaOligosPaginatedList {
  dnaOligos: DnaOligo[];
  nextToken?: string;
}

export interface DnaOligosQueryParams {
  pageSize?: number;
  nextToken?: string;
  sort?: string;
  createdAt?: string;
  modifiedAt?: string;
  name?: string;
  nameIncludes?: string;
  bases?: string;
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
  customNotationId?: string;
}

export interface DnaOligoCreateRequest {
  name: string;
  bases: string;
  aliases?: string[];
  annotations?: DnaAnnotation[];
  authorIds?: string[];
  customFields?: any;
  fields?: any;
  folderId?: string;
  schemaId?: string;
  helm?: string;
}

function randomDnaSequence(length: number): string {
  const dna = 'ATGC';
  let seq = '';
  for (let i = 0; i < length; i++)
    seq += dna[Math.floor(Math.random() * dna.length)];
  return seq;
}

export const mockDnaOligos: DnaOligosPaginatedList = {
  dnaOligos: Array.from({length: 100}, (_, i) => {
    const len = 15 + Math.floor(Math.random() * 16); // 15-30
    return {
      id: `oligo_${String(i+1).padStart(3, '0')}`,
      name: `DNA Oligo ${i+1}`,
      aliases: [`Oligo${i+1}`],
      annotations: [],
      apiURL: `https://benchling.com/api/v2/dna-oligos/oligo_${String(i+1).padStart(3, '0')}`,
      archiveRecord: i % 10 === 0 ? { reason: 'Retired' } : null,
      bases: randomDnaSequence(len),
      createdAt: `2023-${String((i%12)+1).padStart(2, '0')}-01T12:00:00Z`,
      creator: {
        handle: `user${i+1}`,
        id: `ent_${i+1}`,
        name: `User ${i+1}`,
      },
      customFields: {},
      customNotation: null,
      customNotationName: null,
      entityRegistryId: null,
      fields: {},
      folderId: null,
      helm: '',
      length: len,
      modifiedAt: `2023-${String((i%12)+1).padStart(2, '0')}-02T12:00:00Z`,
      nucleotideType: 'DNA',
      registrationOrigin: null,
      registryId: null,
      schema: null,
      webURL: `https://benchling.com/benchling/f/lib_55UxcIps-registry/oligo_${String(i+1).padStart(3, '0')}/edit`,
    };
  }),
  nextToken: undefined,
};

export async function queryDnaOligos(params: DnaOligosQueryParams = {}): Promise<DG.DataFrame> {
  // const token = 'YOUR_BENCHLING_API_TOKEN';
  // const query = buildQueryString(params);
  // const url = `https://benchling.com/api/v2/dna-oligos${query ? '?' + query : ''}`;
  // const response = await grok.dapi.fetchProxy(url, {
  //   headers: {
  //     'Authorization': `Bearer ${token}`,
  //     'Accept': 'application/json',
  //   },
  // });
  // if (!response.ok)
  //   throw new Error(`Benchling API error: ${response.statusText}`);
  // const data = await response.json();
  // const df = dataFrameFromObjects(data.dnaOligos ?? []) ?? DG.DataFrame.create();
  // return df;
  const df = dataFrameFromObjects(mockDnaOligos.dnaOligos);
  return df;
}

export async function postDnaOligo(body: DnaOligoCreateRequest): Promise<DnaOligo> {
  // const token = 'YOUR_BENCHLING_API_TOKEN';
  // const response = await grok.dapi.fetchProxy('https://benchling.com/api/v2/dna-oligos', {
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
  // return data as DnaOligo;
  return mockDnaOligos.dnaOligos[0];
} 
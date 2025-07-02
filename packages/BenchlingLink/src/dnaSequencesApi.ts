// Types for DNA Sequences API
export interface DnaSequence {
  aliases?: string[];
  annotations?: any[];
  apiURL?: string;
  archiveRecord?: any;
  bases: string;
  createdAt?: string;
  creator?: any;
  customFields?: any;
  entityRegistryId?: string | null;
  fields?: any;
  folderId?: string | null;
  id: string;
  isCircular: boolean;
  length: number;
  modifiedAt?: string;
  name: string;
  parts?: any[];
  primers?: any[];
  registrationOrigin?: any;
  registryId?: string | null;
  schema?: any;
  translations?: any[];
  webURL?: string;
}

export const mockDnaSequences: { dnaSequences: DnaSequence[]; nextToken?: string } = {
  dnaSequences: [
    {
      id: 'seq_ABC123',
      name: 'Example DNA 1',
      bases: 'ATGCGTACGTAGCTAGCTAGCTAGCTAGCTA',
      isCircular: false,
      length: 30,
      aliases: ['DNA1'],
      annotations: [],
      apiURL: 'https://benchling.com/api/v2/dna-sequences/seq_ABC123',
      createdAt: '2023-01-01T12:00:00Z',
      modifiedAt: '2023-01-02T12:00:00Z',
      webURL: 'https://benchling.com/benchling/f/lib_55UxcIps-registry/seq_ABC123/edit',
    },
    {
      id: 'seq_DEF456',
      name: 'Example DNA 2',
      bases: 'CGTAGCTAGCTAGCTAGCTAGCTAGCTAGCT',
      isCircular: true,
      length: 31,
      aliases: ['DNA2'],
      annotations: [],
      apiURL: 'https://benchling.com/api/v2/dna-sequences/seq_DEF456',
      createdAt: '2023-02-01T12:00:00Z',
      modifiedAt: '2023-02-02T12:00:00Z',
      webURL: 'https://benchling.com/benchling/f/lib_55UxcIps-registry/seq_DEF456/edit',
    },
  ],
  nextToken: undefined,
};

export const mockDnaSequence: DnaSequence = mockDnaSequences.dnaSequences[0];

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';


export interface DNASequencesQueryParams {
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
}

export async function queryDNASequences(params: DNASequencesQueryParams = {}): Promise<DG.DataFrame> {
  // const token = 'YOUR_BENCHLING_API_TOKEN';
  // const query = buildQueryString(params);
  // const url = `https://benchling.com/api/v2/dna-sequences${query ? '?' + query : ''}`;
  // const response = await grok.dapi.fetchProxy(url, {
  //   headers: {
  //     'Authorization': `Bearer ${token}`,
  //     'Accept': 'application/json',
  //   },
  // });
  // if (!response.ok)
  //   throw new Error(`Benchling API error: ${response.statusText}`);
  // const data = await response.json();
  // const df = DG.DataFrame.fromObjects(data.dnaSequences ?? []) ?? DG.DataFrame.create();
  // return df;
  const df = DG.DataFrame.fromObjects(mockDnaSequences.dnaSequences) ?? DG.DataFrame.create();
  return df;
}

export async function getDNASequenceById(dnaSequenceId: string): Promise<DG.DataFrame> {
  const found = mockDnaSequences.dnaSequences.find((seq) => seq.id === dnaSequenceId);
  const df = DG.DataFrame.fromObjects(found ? [found] : []) ?? DG.DataFrame.create();
  return df;
}

export interface DnaSequenceCreateRequest {
  name: string;
  bases: string;
  aliases?: string[];
  annotations?: any[];
  authorIds?: string[];
  customFields?: any;
  fields?: any;
  folderId?: string;
  schemaId?: string;
}

/**
 * Create a new DNA Sequence via POST /dna-sequences
 * @param body The request body for creating a DNA sequence
 * @returns The created DNA sequence
 */
export async function postDNASequence(body: DnaSequenceCreateRequest): Promise<DnaSequence> {
  // const token = 'YOUR_BENCHLING_API_TOKEN';
  // const response = await grok.dapi.fetchProxy('https://benchling.com/api/v2/dna-sequences', {
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
  return mockDnaSequences.dnaSequences[0] as DnaSequence;
} 
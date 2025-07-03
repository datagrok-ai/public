// Types for DNA Sequences API
import { dataFrameFromObjects } from './utils';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
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
  names_anyOf?: string;
  names_anyOf_caseSensitive?: string;
  creatorIds?: string;
  authorIds_anyOf?: string;
  returning?: string;
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
  // const df = dataFrameFromObjects(data.dnaSequences ?? []) ?? DG.DataFrame.create();
  // return df;
  const df = dataFrameFromObjects(mockDnaSequences.dnaSequences);
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

function randomDnaSequence(length: number): string {
  const dna = 'ATGC';
  let seq = '';
  for (let i = 0; i < length; i++)
    seq += dna[Math.floor(Math.random() * dna.length)];
  return seq;
}

export const mockDnaSequences = {
  dnaSequences: [
    ...Array.from({length: 100}, (_, i) => {
      const len = 20 + Math.floor(Math.random() * 81); // 20-100
      return {
        id: `seq_${String(i+1).padStart(3, '0')}`,
        name: `Example DNA ${i+1}`,
        bases: randomDnaSequence(len),
        isCircular: i % 2 === 0,
        length: len,
        aliases: [`DNA${i+1}`],
        annotations: [],
        apiURL: `https://benchling.com/api/v2/dna-sequences/seq_${String(i+1).padStart(3, '0')}`,
        createdAt: `2023-${String((i%12)+1).padStart(2, '0')}-01T12:00:00Z`,
        modifiedAt: `2023-${String((i%12)+1).padStart(2, '0')}-02T12:00:00Z`,
        webURL: `https://benchling.com/benchling/f/lib_55UxcIps-registry/seq_${String(i+1).padStart(3, '0')}/edit`,
      };
    })
  ],
  nextToken: undefined,
};
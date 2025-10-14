import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { dataFrameFromObjects } from './utils';
import { UserSummary, ArchiveRecord } from './types';

export interface Plate {
  archiveRecord?: ArchiveRecord | null;
  availableCapacity?: number | null;
  barcode?: string | null;
  createdAt?: string;
  creator?: UserSummary;
  fields?: any;
  id: string;
  modifiedAt?: string;
  name: string;
  occupiedCapacity?: number | null;
  parentStorageId?: string | null;
  projectId?: string | null;
  schema?: any;
  totalCapacity?: number | null;
  type: 'matrix_plate' | 'well_plate';
  webURL?: string;
  wells?: Record<string, any>;
}

export interface PlatesPaginatedList {
  plates: Plate[];
  nextToken?: string;
}

export interface PlatesQueryParams {
  pageSize?: number;
  nextToken?: string;
  sort?: string;
  schemaId?: string;
  schemaFields?: string;
  createdAt?: string;
  modifiedAt?: string;
  name?: string;
  nameIncludes?: string;
  emptyPositions?: number;
  emptyPositions_gte?: number;
  emptyPositions_gt?: number;
  emptyPositions_lte?: number;
  emptyPositions_lt?: number;
  emptyContainers?: number;
  emptyContainers_gte?: number;
  emptyContainers_gt?: number;
  emptyContainers_lte?: number;
  emptyContainers_lt?: number;
  ancestorStorageId?: string;
  storageContentsId?: string;
  storageContentsIds?: string;
  archiveReason?: string;
  ids?: string;
  barcodes?: string;
  names_anyOf?: string;
  names_anyOf_caseSensitive?: string;
  returning?: string;
  creatorIds?: string;
}

export interface PlateCreateRequest {
  barcode?: string;
  containerSchemaId?: string;
  fields?: any;
  name: string;
  parentStorageId?: string;
  projectId?: string;
  schemaId: string;
  wells?: Record<string, any>;
}

export const mockPlates: PlatesPaginatedList = {
  plates: Array.from({length: 100}, (_, i) => ({
    id: `plt_${String(i+1).padStart(3, '0')}`,
    name: `Plate ${i+1}`,
    barcode: `W${100000 + i}`,
    createdAt: `2023-${String((i%12)+1).padStart(2, '0')}-01T12:00:00Z`,
    modifiedAt: `2023-${String((i%12)+1).padStart(2, '0')}-02T12:00:00Z`,
    type: i % 2 === 0 ? 'matrix_plate' : 'well_plate',
    availableCapacity: 96,
    occupiedCapacity: Math.floor(Math.random() * 96),
    totalCapacity: 96,
    parentStorageId: null,
    projectId: null,
    archiveRecord: i % 10 === 0 ? { reason: 'Retired' } : null,
    creator: {
      handle: `user${i+1}`,
      id: `ent_${i+1}`,
      name: `User ${i+1}`,
    },
    webURL: `https://benchling.com/plates/plt_${String(i+1).padStart(3, '0')}`,
    wells: {},
    fields: {},
    schema: null,
  })),
  nextToken: undefined,
};

export async function queryPlates(params: PlatesQueryParams = {}): Promise<DG.DataFrame> {
  // const token = 'YOUR_BENCHLING_API_TOKEN';
  // const query = buildQueryString(params);
  // const url = `https://benchling.com/api/v2/plates${query ? '?' + query : ''}`;
  // const response = await grok.dapi.fetchProxy(url, {
  //   headers: {
  //     'Authorization': `Bearer ${token}`,
  //     'Accept': 'application/json',
  //   },
  // });
  // if (!response.ok)
  //   throw new Error(`Benchling API error: ${response.statusText}`);
  // const data = await response.json();
  // const df = dataFrameFromObjects(data.plates ?? []) ?? DG.DataFrame.create();
  // return df;
  const df = dataFrameFromObjects(mockPlates.plates);
  return df;
}

export async function postPlate(body: PlateCreateRequest): Promise<Plate> {
  // const token = 'YOUR_BENCHLING_API_TOKEN';
  // const response = await grok.dapi.fetchProxy('https://benchling.com/api/v2/plates', {
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
  // return data as Plate;
  return mockPlates.plates[0];
} 
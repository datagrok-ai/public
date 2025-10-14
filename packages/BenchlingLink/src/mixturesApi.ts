import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { dataFrameFromObjects } from './utils';
import { UserSummary, ArchiveRecord } from './types';

export interface Ingredient {
  amount?: string | null;
  catalogIdentifier?: string | null;
  componentEntity?: any;
  componentLotContainer?: any;
  componentLotEntity?: any;
  componentLotText?: string | null;
  hasParent?: boolean;
  notes?: string | null;
  targetAmount?: string | null;
  units?: string | null;
}

export interface Mixture {
  aliases?: string[];
  allowMeasuredIngredients?: boolean;
  amount?: string;
  apiURL?: string;
  archiveRecord?: ArchiveRecord | null;
  authors?: UserSummary[];
  createdAt?: string;
  creator?: UserSummary;
  customFields?: any;
  entityRegistryId?: string | null;
  fields?: any;
  folderId?: string | null;
  id: string;
  ingredients?: Ingredient[];
  modifiedAt?: string;
  name: string;
  registrationOrigin?: any;
  registryId?: string | null;
  schema?: any;
  units?: string | null;
  webURL?: string;
}

export interface MixturesPaginatedList {
  mixtures: Mixture[];
  nextToken?: string;
}

export interface MixturesQueryParams {
  nextToken?: string;
  pageSize?: number;
  sort?: string;
  createdAt?: string;
  modifiedAt?: string;
  name?: string;
  nameIncludes?: string;
  folderId?: string;
  mentionedIn?: string;
  projectId?: string;
  registryId?: string;
  schemaId?: string;
  schemaFields?: string;
  archiveReason?: string;
  mentions?: string;
  ids?: string;
  names_anyOf?: string;
  names_anyOf_caseSensitive?: string;
  entityRegistryIds_anyOf?: string;
  ingredientComponentEntityIds?: string;
  ingredientComponentEntityIds_anyOf?: string;
  authorIds_anyOf?: string;
}

export interface MixtureCreateRequest {
  name: string;
  ingredients: Ingredient[];
  schemaId: string;
  units: string;
  aliases?: string[];
  amount?: string;
  authorIds?: string[];
  customFields?: any;
  entityRegistryId?: string;
  fields?: any;
  folderId?: string;
}

export const mockMixtures: MixturesPaginatedList = {
  mixtures: Array.from({length: 100}, (_, i) => ({
    id: `mxt_${String(i+1).padStart(3, '0')}`,
    name: `Mixture ${i+1}`,
    aliases: [`Mix${i+1}`],
    allowMeasuredIngredients: true,
    amount: `${100 + i}`,
    apiURL: `https://benchling.com/api/v2/mixtures/mxt_${String(i+1).padStart(3, '0')}`,
    archiveRecord: i % 10 === 0 ? { reason: 'Retired' } : null,
    authors: [
      {
        handle: `user${i+1}`,
        id: `ent_${i+1}`,
        name: `User ${i+1}`,
      }
    ],
    createdAt: `2023-${String((i%12)+1).padStart(2, '0')}-01T12:00:00Z`,
    creator: {
      handle: `user${i+1}`,
      id: `ent_${i+1}`,
      name: `User ${i+1}`,
    },
    customFields: {},
    entityRegistryId: null,
    fields: {},
    folderId: null,
    ingredients: [
      {
        amount: `${10 + i}`,
        catalogIdentifier: `CAT${i+1}`,
        componentEntity: { id: `ent_${i+1}` },
        componentLotContainer: null,
        componentLotEntity: null,
        componentLotText: null,
        hasParent: false,
        notes: `Note for ingredient ${i+1}`,
        targetAmount: null,
        units: 'mg',
      }
    ],
    modifiedAt: `2023-${String((i%12)+1).padStart(2, '0')}-02T12:00:00Z`,
    registrationOrigin: null,
    registryId: null,
    schema: null,
    units: 'mg',
    webURL: `https://benchling.com/benchling/f/lib_55UxcIps-registry/mxt_${String(i+1).padStart(3, '0')}/edit`,
  })),
  nextToken: undefined,
};

export async function queryMixtures(params: MixturesQueryParams = {}): Promise<DG.DataFrame> {
  // const token = 'YOUR_BENCHLING_API_TOKEN';
  // const query = buildQueryString(params);
  // const url = `https://benchling.com/api/v2/mixtures${query ? '?' + query : ''}`;
  // const response = await grok.dapi.fetchProxy(url, {
  //   headers: {
  //     'Authorization': `Bearer ${token}`,
  //     'Accept': 'application/json',
  //   },
  // });
  // if (!response.ok)
  //   throw new Error(`Benchling API error: ${response.statusText}`);
  // const data = await response.json();
  // const df = dataFrameFromObjects(data.mixtures ?? []) ?? DG.DataFrame.create();
  // return df;
  const df = dataFrameFromObjects(mockMixtures.mixtures);
  return df;
}

export async function postMixture(body: MixtureCreateRequest): Promise<Mixture> {
  // const token = 'YOUR_BENCHLING_API_TOKEN';
  // const response = await grok.dapi.fetchProxy('https://benchling.com/api/v2/mixtures', {
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
  // return data as Mixture;
  return mockMixtures.mixtures[0];
} 
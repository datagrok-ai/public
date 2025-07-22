// Types for Molecules API
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { dataFrameFromObjects } from './utils';
import { _package } from './package';
import { UserSummary } from './types';

export interface Molecule {
  aliases?: string[];
  apiURL?: string;
  archiveRecord?: any;
  canonicalizedSmiles?: string;
  createdAt?: string;
  creator?: UserSummary;
  customFields?: any;
  entityRegistryId?: string | null;
  fields?: any;
  folderId?: string | null;
  id: string;
  modifiedAt?: string;
  name: string;
  originalSmiles?: string | null;
  registrationOrigin?: any;
  registryId?: string | null;
  schema?: any;
  smiles: string;
  webURL?: string;
}

export interface MoleculesPaginatedList {
  molecules: Molecule[];
  nextToken?: string;
}

export async function getSmilesListFromCsv(): Promise<DG.Column> {
  const df = await DG.DataFrame.fromCsv(await _package.files.readAsText('benchling_mock_molecules.csv'));
  return df.col('canonical_smiles')!;
}

export async function getMockMolecules(): Promise<MoleculesPaginatedList> {
  const smilesList = await getSmilesListFromCsv();
  return {
    molecules: Array.from({length: 100}, (_, i) => {
      const smiles = smilesList.get(i);
      return {
        id: `mol_${String(i+1).padStart(3, '0')}`,
        name: `Example Molecule ${i+1}`,
        aliases: [`Mol${i+1}`],
        apiURL: `https://benchling.com/api/v2/molecules/mol_${String(i+1).padStart(3, '0')}`,
        archiveRecord: null,
        canonicalizedSmiles: smiles,
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
        modifiedAt: `2023-${String((i%12)+1).padStart(2, '0')}-02T12:00:00Z`,
        originalSmiles: smiles,
        registrationOrigin: null,
        registryId: null,
        schema: null,
        smiles: smiles,
        webURL: `https://benchling.com/benchling/f/lib_55UxcIps-registry/mol_${String(i+1).padStart(3, '0')}/edit`,
      };
    }),
    nextToken: undefined,
  };
}

export async function getMockMolecule(): Promise<Molecule> {
  return (await getMockMolecules()).molecules[0];
}

export interface MoleculesQueryParams {
  pageSize?: number;
  nextToken?: string;
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
  entityRegistryIds_anyOf?: string;
  names_anyOf?: string;
  authorIds_anyOf?: string;
  chemicalSubstructure_mol?: string;
  chemicalSubstructure_smiles?: string;
}

export async function queryMolecules(params: MoleculesQueryParams = {}): Promise<DG.DataFrame> {
  // const token = 'YOUR_BENCHLING_API_TOKEN';
  // const query = buildQueryString(params);
  // const url = `https://benchling.com/api/v2/molecules${query ? '?' + query : ''}`;
  // const response = await grok.dapi.fetchProxy(url, {
  //   headers: {
  //     'Authorization': `Bearer ${token}`,
  //     'Accept': 'application/json',
  //   },
  // });
  // if (!response.ok)
  //   throw new Error(`Benchling API error: ${response.statusText}`);
  // const data = await response.json();
  // const df = dataFrameFromObjects(data.molecules ?? []) ?? DG.DataFrame.create();
  // return df;
  const df = dataFrameFromObjects(await getMockMolecules().then(m => m.molecules));
  return df;
}

export interface MoleculeCreateRequest {
  name: string;
  smiles: string;
  formula?: string;
}

export async function postMolecule(body: MoleculeCreateRequest): Promise<Molecule> {
  // const token = 'YOUR_BENCHLING_API_TOKEN';
  // const response = await grok.dapi.fetchProxy('https://benchling.com/api/v2/molecules', {
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
  // return data as Molecule;
  return (await getMockMolecules()).molecules[0];
} 
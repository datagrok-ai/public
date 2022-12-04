import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getRdKitModule, getRdKitWebRoot} from '../utils/chem-common-rdkit';
import {_convertMolNotation} from '../utils/convert-notation-utils';

let unichemSources: DG.DataFrame;

export async function getOrLoadUnichemSources(): Promise<DG.DataFrame> {
  unichemSources ??= await grok.data.loadTable(`${getRdKitWebRoot()}files/unichem-sources.csv`);
  return unichemSources;
}

class UniChemSource {
  id: number;
  name: string;
  fullName: string;
  labelName: string;
  baseUrl: string;
  homePage: string;
  description: string;

  static _sources: {[key: number]: UniChemSource} = {};

  static idNames : {[key: number]: string} = {
    1: 'chembl',
    2: 'drugbank',
    3: 'pdb',
    4: 'gtopdb',
    5: 'pubchem_dotf',
    6: 'kegg_ligand',
    7: 'chebi',
    8: 'nih_ncc',
    9: 'zinc',
    10: 'emolecules',
    11: 'ibm',
    12: 'atlas',
    14: 'fdasrs',
    15: 'surechembl',
    17: 'pharmgkb',
    18: 'hmdb',
    20: 'selleck',
    21: 'pubchem_tpharma',
    22: 'pubchem',
    23: 'mcule',
    24: 'nmrshiftdb2',
    25: 'lincs',
    26: 'actor',
    27: 'recon',
    28: 'molport',
    29: 'nikkaji',
    31: 'bindingdb',
    32: 'comptox',
    33: 'lipidmaps',
    34: 'drugcentral',
    35: 'carotenoiddb',
    36: 'metabolights',
    37: 'brenda',
    38: 'rhea',
    39: 'chemicalbook',
    41: 'swisslipids',
    43: 'gsrs',
  };

  constructor(
    id: number, name: string, fullName:string, labelName: string, baseUrl: string, homePage: string,
    description: string,
  ) {
    this.id = id;
    this.name = name;
    this.fullName = fullName;
    this.labelName = labelName;
    this.baseUrl = baseUrl;
    this.homePage = homePage;
    this.description = description;
  }

  static async refreshSources(): Promise<void> {
    const table = await getOrLoadUnichemSources();
    const rowCount = table.rowCount;
    for (let i = 0; i < rowCount; i++) {
      const id = table.get('src_id', i);
      this._sources[id] = new UniChemSource(
        id, table.get('name', i), table.get('name_long', i), table.get('name_label', i),
        table.get('base_id_url', i), table.get('base_id_url', i), table.get('description', i),
      );
    }
  }

  static byName(name: string): UniChemSource | null {
    for (const source of Object.values(this._sources)) {
      if (source.name === name)
        return source;
    }
    return null;
  }
}

async function getCompoundsIds(inchiKey: string): Promise<{[k:string]: any} | undefined> {
  const url = `https://www.ebi.ac.uk/unichem/rest/inchikey/${inchiKey}`;
  const params: RequestInit = {method: 'GET', referrerPolicy: 'strict-origin-when-cross-origin'};
  const response = await grok.dapi.fetchProxy(url, params);
  const json = await response.json();
  if (json.error)
    return;
  const sources: {[key: string]: any}[] = json.filter((s: {[key: string]: string | number}) => {
    const srcId = parseInt(`${s['src_id']}`);
    s['src_id'] = srcId;
    return srcId in UniChemSource.idNames;
  });

  return response.status !== 200 ?
    {} : Object.fromEntries(sources.map((m) => [UniChemSource.idNames[(m['src_id'] as number)], m['src_compound_id']]));
}

export async function getIdMap(smiles: string): Promise<{[k:string]: any} | null> {
  const rdKitModule = getRdKitModule();
  const mol = rdKitModule.get_mol(smiles);
  const inchiKey = rdKitModule.get_inchikey_for_inchi(mol.get_inchi());
  mol.delete();
  const idMap = await getCompoundsIds(inchiKey);

  if (typeof idMap === 'undefined')
    return null;

  await UniChemSource.refreshSources();

  for (const [source, id] of Object.entries(idMap))
    idMap[source] = {id: id, link: UniChemSource.byName(source)!.baseUrl + id};
    // idMap[source] = ui.link(id, () => window.open(UniChemSource.byName(source)!.baseUrl + id));

  return idMap;
}

export async function identifiersWidget(smiles: string): Promise<DG.Widget> {
  const rdKitModule = getRdKitModule();
  try {
    smiles = _convertMolNotation(smiles, 'unknown', 'smiles', rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possible malformed'));
  }
  let idMap: {[k: string]: any} | null = null;
  try {
    idMap = await getIdMap(smiles);
  } catch (e) {
    console.warn(e);
  }
  if (idMap === null)
    return new DG.Widget(ui.divText('Not found in UniChem'));
  for (const [source, identifier] of Object.entries(idMap))
    idMap[source] = ui.link(identifier.id, () => window.open(identifier.link));
  return new DG.Widget(ui.tableFromMap(idMap));
}

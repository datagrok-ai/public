import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getRdKitModule, getRdKitWebRoot} from '../utils/chem-common-rdkit';
import {_convertMolNotation} from '../utils/convert-notation-utils';

const CHEMBL = 'Chembl';
const PUBCHEM = 'PubChem';

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
    1: CHEMBL.toLowerCase(),
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
    22: PUBCHEM.toLowerCase(),
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

export async function getIdMap(inchiKey: string): Promise<{[k:string]: any} | null> {
  const idMap = await getCompoundsIds(inchiKey);

  if (typeof idMap === 'undefined')
    return null;

  await UniChemSource.refreshSources();

  for (const [source, id] of Object.entries(idMap))
    idMap[source] = {id: id, link: UniChemSource.byName(source)!.baseUrl + id};
    // idMap[source] = ui.link(id, () => window.open(UniChemSource.byName(source)!.baseUrl + id));

  return idMap;
}

async function getIUPACName(smiles: string): Promise<string> {
  const preparedSmiles = smiles.replaceAll(`#`, `%23`); // need to escape # sign (triple bond) in URL
  const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${preparedSmiles}/property/IUPACName/JSON`;
  const response = await fetch(url);
  const responseJson = await response.json();
  const result = responseJson.PropertyTable?.Properties;
  return (result && result[0].hasOwnProperty('IUPACName')) ? result[0].IUPACName : 'Not found in PubChem';
}

export async function identifiers(molfile: string): Promise<HTMLElement> {
  const rdKitModule = getRdKitModule();
  const mol = rdKitModule.get_mol(molfile);
  const inchi = mol.get_inchi();
  const inchiKey = rdKitModule.get_inchikey_for_inchi(inchi);
  mol.delete();

  let idMap: {[k: string]: any} | null = null;
  try {
    idMap = await getIdMap(inchiKey);
  } catch (e) {
    console.warn(e);
  }

  const identifiersDiv = ui.div();
  const mainIdentifiers = createMainIdentifiersMap(molfile, inchi, inchiKey, idMap);
  identifiersDiv.append(ui.tableFromMap(mainIdentifiers));
  
  if (idMap === null)
    identifiersDiv.append(ui.divText('Not found in UniChem'));
  else {
    const otherIdentifiers: {[k: string]: any} = {};
    for (const [source, identifier] of Object.entries(idMap))
    otherIdentifiers[source] = ui.link(identifier.id, () => window.open(identifier.link));
    const acc = ui.accordion();
    acc.addPane(`Other identifiers`, 
      () => Object.keys(otherIdentifiers).length === 0 ? 
        ui.divText('No otehr identifiers found') :  ui.tableFromMap(otherIdentifiers));
    identifiersDiv.append(acc.root);
  }

  return identifiersDiv;
}


function createMainIdentifiersMap(molfile: string, inchi: string, inchiKey: string, idMap: {[k: string]: any} | null): {[k: string]: any} {
  const map: {[k: string]: any} = {};
  const smiles =  _convertMolNotation(molfile, DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles, getRdKitModule());
  map['Name'] = ui.wait(async () => ui.divText(await getIUPACName(smiles)));
  map['Smiles'] = smiles;
  map['Inchi'] = inchi;
  map['Inchi key'] = inchiKey;
  extractMainIdentifier(CHEMBL, map, idMap);
  extractMainIdentifier(PUBCHEM, map, idMap);
  function extractMainIdentifier(source: string, map: {[k: string]: any}, idMap: {[k: string]: any} | null) {
    if (idMap) {
      const identifier = idMap[source.toLowerCase()];
      if (identifier) {
        map[source] = ui.link(identifier.id, () => window.open(identifier.link));
        delete idMap[source.toLowerCase()];
      } else {
        map[source] = '-';
      }
    }
  }
  return map;
}
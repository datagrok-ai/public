import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getRdKitModule, getRdKitWebRoot} from '../utils/chem-common-rdkit';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {getMolSafe} from '../utils/mol-creation_rdkit';
import { checkPackage } from '../utils/elemental-analysis-utils';
import {inchiToInchiKey, inchiToSmiles, smilesToCanonical, smilesToInchi, smilesToInchiKey} from "../docker/api";

const CHEMBL = 'Chembl';
const PUBCHEM = 'PubChem';
const MW_DESC = 'exactmw';
const MAX_MW_FOR_IDENTIFIERS = 1200;

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

  static idNamesChoices: string = `actor,atlas,bindingdb,brenda,carotenoiddb,chebi,chembl,chemicalbook,comptox,drugbank,drugcentral,emolecules,fdasrs,gtopdb,hmdb,ibm,inchi,inchi_key,kegg_ligand,lincs,lipidmaps,mcule,metabolights,molport,nih_ncc,nikkaji,nmrshiftdb2,pdb,pharmgkb,pubchem,pubchem_dotf,pubchem_tpharma,recon,rhea,selleck,smiles,surechembl,zinc`;

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

async function getCompoundsIds(inchiKey: string): Promise<{[k:string]: any}> {
  const json = await grok.functions.call(`ChemblApi:getCompoundsIds`, { 'inchiKey': inchiKey });
  const sources: {[key: string]: any}[] = json.filter((s: {[key: string]: string | number}) => {
    const srcId = parseInt(`${s['src_id']}`);
    s['src_id'] = srcId;
    return srcId in UniChemSource.idNames;
  });
  return Object.fromEntries(sources.map((m) => [UniChemSource.idNames[(m['src_id'] as number)], m['src_compound_id']]));
}

export async function getIdMap(inchiKey: string): Promise<{[k:string]: any} | null> {
  const idMap = await getCompoundsIds(inchiKey);

  await UniChemSource.refreshSources();

  for (const [source, id] of Object.entries(idMap))
    idMap[source] = {id: id, link: UniChemSource.byName(source)!.baseUrl + id};
    // idMap[source] = ui.link(id, () => window.open(UniChemSource.byName(source)!.baseUrl + id));

  return idMap;
}

export async function identifiersWidget(molfile: string): Promise<DG.Widget> {
  const rdKitModule = getRdKitModule();
  const mol = getMolSafe(molfile, {}, rdKitModule);
  let fullIdentifiersPanel = true;
  if (mol.mol) {
    const inchi = mol.mol.get_inchi();
    const inchiKey = rdKitModule.get_inchikey_for_inchi(inchi);
    if (JSON.parse(mol.mol.get_descriptors())[MW_DESC] > MAX_MW_FOR_IDENTIFIERS)
      fullIdentifiersPanel = false;
    mol.mol.delete();

    let idMap: {[k: string]: any} | null = null;
    try {
      if (checkPackage('ChemblApi', 'getCompoundsIds'))
        idMap = await getIdMap(inchiKey);
    } catch (e) {
      console.warn(e);
    }

    const mainIdentifiers = createIdentifiersMap(molfile, inchi, inchiKey, idMap, fullIdentifiersPanel);
    return new DG.Widget(ui.tableFromMap(mainIdentifiers));
  } else
    return new DG.Widget(ui.divText('Malformed molecule'));
}


function createIdentifiersMap(molfile: string, inchi: string, inchiKey: string,
  idMap: {[k: string]: any} | null, fullIdentifiersPanel: boolean): {[k: string]: any} {
  const map: {[k: string]: any} = {};
  const smiles = _convertMolNotation(molfile, DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles, getRdKitModule());
  if (fullIdentifiersPanel)
    addNameField(map, smiles);
  map['Smiles'] = smiles;
  map['Inchi'] = inchi;
  map['Inchi key'] = inchiKey;
  if (fullIdentifiersPanel) {
    if (idMap && Object.keys(idMap).length) {
      extractMainIdentifier(CHEMBL, map, idMap);
      extractMainIdentifier(PUBCHEM, map, idMap);
      for (const [source, identifier] of Object.entries(idMap))
        map[source] = ui.link(identifier.id, () => window.open(identifier.link));
    }
    function extractMainIdentifier(source: string, map: {[k: string]: any}, idMap: {[k: string]: any} | null) {
      const identifier = idMap ? idMap[source.toLowerCase()] : null;
      if (identifier) {
        map[source] = ui.link(identifier.id, () => window.open(identifier.link));
        delete idMap![source.toLowerCase()];
      } else
        map[source] = '-';
    }
  }
  return map;
}

function addNameField(map: { [key: string]: any }, smiles: string) {
  const packageExists = checkPackage('PubchemApi', 'GetIupacName');
  if (packageExists)
    map['Name'] = ui.wait(async () => ui.divText(await grok.functions.call(`PubChemApi:GetIupacName`, { 'smiles': smiles })));
}

export async function openMapIdentifiersDialog() {
  const dlg = ui.dialog('Map identifiers...');
  if (!grok.shell.t) {
    grok.shell.warning('There are no table');
    return;
  }
  const table = grok.shell.t;
  const columnSelector = ui.columnInput('Molecules', table,
    table.columns.bySemTypeAll(DG.SEMTYPE.MOLECULE).find((_) => true) ?? null, null,
    {filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MOLECULE});
  columnSelector.root.children[0].classList.add('d4-chem-descriptors-molecule-column-input');
  const choices = UniChemSource.idNamesChoices.split(',');
  const fromSource = ui.choiceInput('fromSource', 'smiles', choices, null, {nullable: false});
  const toSource = ui.choiceInput('toSource', choices[0], choices, null, {nullable: false});
  dlg.add(columnSelector);
  dlg.add(fromSource);
  dlg.add(toSource);
  dlg.onOK(async () => {
    if (!columnSelector.value) {
      grok.shell.warning('Column should be chosen');
      return;
    }
    await getMapIdentifiers(table, columnSelector.value, fromSource.value!, toSource.value!);
  });
  dlg.show();
}

export async function getMapIdentifiers(table: DG.DataFrame, ids: DG.Column, fromSource: string,
                                               toSource: string) {
  let result: string[] | DG.Column | null;
  const molecules = ids.toList();
  if (fromSource === 'smiles') {
    if (toSource === 'smiles')
      result = await smilesToCanonical(molecules);
    else if (toSource === 'inchi')
      result = await smilesToInchi(molecules);
    else if (toSource === 'inchi_key')
      result = await smilesToInchiKey(molecules);
    else {
      let d = await _chemMapViaQuery(DG.Column.fromList('string', 'key',
        await smilesToInchiKey(molecules)), 'inchiKeyToChembl', 'chembl');
      if (toSource !== 'chembl') {
        d = await _chemMapIdentifiersUnichem(d, 'chembl', toSource);
        result = d.getCol(toSource);
      }
      else
        result = d.getCol('chembl');
      d.columns.remove(result.name);
    }
  }
  else if (fromSource === 'inchi') {
    if (toSource === 'smiles')
      result = await inchiToSmiles(molecules);
    else if (toSource === 'inchi_key')
      result = await inchiToInchiKey(molecules);
    else {
      let d = await _chemMapViaQuery(DG.Column.fromList('string', 'key', await inchiToInchiKey(molecules)), 'inchiKeyToChembl', 'chembl');
      if (toSource !== 'chembl') {
        d = await _chemMapIdentifiersUnichem(d, 'chembl', toSource);
        result = d.getCol(toSource);
      }
      else
        result = d.getCol('chembl');
    }
  }
  else if (fromSource == 'inchi_key') {
    if (toSource === 'smiles')
      result = (await _chemMapViaQuery(ids, 'inchiKeyToSmiles', toSource))?.col(toSource);
    else if (toSource === 'inchi')
      result = (await _chemMapViaQuery(ids, 'inchiKeyToInchi', toSource))?.col(toSource);
    else {
      let d = await _chemMapViaQuery(ids, 'inchiKeyToChembl', 'chembl');
      if (toSource !== 'chembl') {
        d = await _chemMapIdentifiersUnichem(d, 'chembl', toSource);
        result = d.getCol(toSource);
      }
      else
        result = d.getCol('chembl');
    }
  }
  else {
    async function unichemToSourceViaQuery(queryName: string)  {
      let res = await _chemMapIdentifiersUnichem(table, fromSource, 'chembl');
      if (res.rowCount !== 0)
        result = (await _chemMapViaQuery(res.getCol(res.columns.names()[0]), queryName, toSource)).getCol(toSource);
    }
    if (toSource == 'smiles')
      return await unichemToSourceViaQuery('chemblToSmiles');
    else if (toSource == 'inchi')
      return await unichemToSourceViaQuery('chemblToInchi');
    else if (toSource == 'inchi_key')
      return await unichemToSourceViaQuery('chemblToInchiKey');
    else
      return await _chemMapIdentifiersUnichem(table, fromSource, toSource);
  }
  if (result instanceof DG.Column) {
    result.dataFrame.columns.remove(result.name);
    result.name = table.columns.getUnusedName(result.name);
    result.setTag('MoleculeId', toSource);
    table.columns.add(result);
  }
}

async function _chemMapViaQuery(keys: DG.Column, queryName: string, resultColumnName: string) {
  const query: DG.DataQuery = await grok.functions.eval(`chembl:${queryName}`);
  if (!query) throw new Error('Missing Chembl package');
  const df = DG.DataFrame.fromColumns([keys]);
  df.columns.add(DG.Column.fromList('int', 'order', [...Array(df.rowCount).keys()]))
  const call: DG.FuncCall = await (query.prepare({'ids' : df})).call();
  const result: DG.DataFrame = call.getOutputParamValue() as DG.DataFrame;
  if (result.columns.length === 0) return result;
  result.getCol(result.columns.names()[0]).name = resultColumnName;
  return result;
}

async function _chemMapIdentifiersUnichem(table: DG.DataFrame, fromSource: string, toSource: string): Promise<DG.DataFrame> {
  if (fromSource === toSource)
    return table;
  let fromId = UniChemSource.byName(fromSource)!.id;
  let toId = UniChemSource.byName(toSource)!.id;
  const conn: DG.DataConnection = await grok.functions.eval('chembl:unichem');
  table.getCol(table.columns.names()[0]).name = 'keys'
  let n: number = fromId;
  let m: number = toId;
  if (n > m) {
    n = toId;
    m = fromId;
  }
  let tableName = `src${n}src${m}`;
  table.columns.add(DG.Column.fromList('int', 'order', [...Array(table.rowCount).keys()]));
  const dataQuery: DG.DataQuery = conn.query('', `--input: dataframe ids\n
select (select ${n == fromId ? 'to_id' : 'from_id'} from ${tableName}
where ids.keys = ${tableName}.${n == fromId ? 'from_id' : 'to_id'} limit 1) from  
ids order by ids.order`);
  const call = await (dataQuery.prepare({'ids': table})).call()
  const result: DG.DataFrame = call.getOutputParamValue() as DG.DataFrame;
  if (result.columns.length === 0) return result;
  result.getCol(result.columns.names()[0]).name = toSource;
  return result;
}

/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getBy} from './pubchem';
import {COLUMN_NAMES, getSearchWidget} from './widget';
import {pubChemRest} from './tests/const';
import {TRIPLE_BOND, TRIPLE_BOND_REPLACE_SYMBOL} from './constants';
export * from './package.g';
export const _package = new DG.Package();

/*
name: PubChem | Info
tags: panel, widgets
input: string molString {semType: Molecule}
output: widget result
export async function pubChemPanel(molString: string): Promise<DG.Widget> {
  const pubChemId = await smilesToPubChem(molString);
  return new DG.Widget(ui.wait(async () => await buildAccordion(pubChemId)));
}
*/

export class PackageFunctions {
  @grok.decorators.panel({
    'name': 'Databases | PubChem | Substructure Search',
    'meta': {role: 'widgets'},
  })
  static async pubChemSubstructureSearchPanel(
    @grok.decorators.param({'options': {'semType': 'Molecule'}}) molString: string): Promise<DG.Widget> {
    return molString ? await getSearchWidget(molString, 'substructure') : new DG.Widget(ui.divText('SMILES is empty'));
  }


  @grok.decorators.panel({
    'name': 'Databases | PubChem | Similarity Search',
    'meta': {role: 'widgets'},
  })
  static async pubChemSimilaritySearchPanel(
    @grok.decorators.param({'options': {'semType': 'Molecule'}}) molString: string): Promise<DG.Widget> {
    return molString ? await getSearchWidget(molString, 'similarity') : new DG.Widget(ui.divText('SMILES is empty'));
  }


  @grok.decorators.panel({
    'name': 'Databases | PubChem | Identity Search',
    'meta': {role: 'widgets'},
  })
  static async pubChemIdentitySearch(
    @grok.decorators.param({'options': {'semType': 'Molecule'}}) molString: string): Promise<DG.Widget> {
    return molString ? await getSearchWidget(molString, 'identity') : new DG.Widget(ui.divText('SMILES is empty'));
  }


  @grok.decorators.func({
    'meta': {
      'role': 'converter',
      'inputRegexp': '(^\\s*[Pp][Uu][Bb][Cc][Hh][Ee][Mm]\\s*\\:\\s*[0-9]+\\s*$)',
    },
    'connection': 'PubChemApi',
    'outputs': [{name: 'result', type: 'string', options: {semType: 'Molecule'}}],
  })
  static async pubChemToSmiles(
    id: string): Promise<string> {
    const pubChemId = id.substring(id.indexOf(':') + 1).trim();
    const url = `${pubChemRest}/pug/compound/cid/${pubChemId}/property/CanonicalSMILES/JSON`;
    const response = await grok.dapi.fetchProxy(url);
    const json = await response.json();
    return json['PropertyTable']['Properties'][0][COLUMN_NAMES.CANONICAL_SMILES] ??
      json['PropertyTable']['Properties'][0][COLUMN_NAMES.CONNECTIVITY_SMILES];
  }


  @grok.decorators.func({
    'meta': {
      'role': 'converter',
      'inputRegexp': '([A-Z]{14}-[A-Z]{10}-N)',
    },
    'connection': 'PubChemApi',
    'outputs': [{name: 'result', type: 'string', options: {semType: 'Molecule'}}],
  })
  static async inchiKeysToSmiles(
    id: string): Promise<string> {
    const s = await getBy('InChIKey', 'cids', id);
    const cids = s['IdentifierList']['CID'][0];
    const smiles = await PackageFunctions.pubChemToSmiles(cids.toString());
    return smiles;
  }


  @grok.decorators.func({
    'connection': 'PubChemApi',
  })
  static async GetIupacName(
    smiles: string): Promise<string> {
    // need to escape # sign (triple bond) in URL
    const preparedSmiles = smiles.replaceAll(TRIPLE_BOND, TRIPLE_BOND_REPLACE_SYMBOL);
    const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${preparedSmiles}/property/IUPACName/JSON`;
    const response = await grok.dapi.fetchProxy(url);
    const responseJson = await response.json();
    const result = responseJson.PropertyTable?.Properties;
    return (result && result[0].hasOwnProperty('IUPACName')) ? result[0].IUPACName : 'Not found in PubChem';
  }
}

/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getBy} from './pubchem';
import {getSearchWidget} from './widget';
import {pubChemRest} from './tests/const';
import { TRIPLE_BOND, TRIPLE_BOND_REPLACE_SYMBOL } from './constants';

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

//name: Databases | PubChem | Substructure Search
//tags: panel, widgets
//input: string molString {semType: Molecule}
//output: widget result
export async function pubChemSubstructureSearchPanel(molString: string): Promise<DG.Widget> {
  return molString ? await getSearchWidget(molString, 'substructure') : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Databases | PubChem | Similarity Search
//tags: panel, widgets
//input: string molString {semType: Molecule}
//output: widget result
export async function pubChemSimilaritySearchPanel(molString: string): Promise<DG.Widget> {
  return molString ? await getSearchWidget(molString, 'similarity') : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Databases | PubChem | Identity Search
//tags: panel, widgets
//input: string molString {semType: Molecule}
//output: widget result
export async function pubChemIdentitySearch(molString: string): Promise<DG.Widget> {
  return molString ? await getSearchWidget(molString, 'identity') : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: pubChemToSmiles
//input: string id
//output: string smiles {semType: Molecule}
//meta.role: converter
//meta.inputRegexp: (^\s*[Pp][Uu][Bb][Cc][Hh][Ee][Mm]\s*\:\s*[0-9]+\s*$)
//connection: PubChemApi
export async function pubChemToSmiles(id: string) {
  const pubChemId = id.substring(id.indexOf(":") + 1).trim();
  const url = `${pubChemRest}/pug/compound/cid/${pubChemId}/property/CanonicalSMILES/JSON`;
  const response = await grok.dapi.fetchProxy(url);
  const json = await response.json();
  return json['PropertyTable']['Properties'][0]['CanonicalSMILES'];
}

//name: inchiKeysToSmiles
//input: string id
//output: string smiles {semType: Molecule}
//meta.role: converter
//meta.inputRegexp: ([A-Z]{14}-[A-Z]{10}-N)
//connection: PubChemApi
export async function inchiKeysToSmiles(id: string) {
  const s = await getBy('InChIKey', 'cids', id);
  const cids = s['IdentifierList']['CID'][0];
  const smiles = await pubChemToSmiles(cids.toString());
  return smiles;
}

//name: GetIupacName
//input: string smiles
//output: string name
//connection: PubChemApi
export async function GetIupacName(smiles: string) {
    // need to escape # sign (triple bond) in URL
    const preparedSmiles = smiles.replaceAll(TRIPLE_BOND, TRIPLE_BOND_REPLACE_SYMBOL);
    const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${preparedSmiles}/property/IUPACName/JSON`;
    const response = await grok.dapi.fetchProxy(url);
    const responseJson = await response.json();
    const result = responseJson.PropertyTable?.Properties;
    return (result && result[0].hasOwnProperty('IUPACName')) ? result[0].IUPACName : 'Not found in PubChem';
}

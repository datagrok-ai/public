/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {init, smilesToPubChem} from './pubchem';
import {pubChemSearchWidget} from './widget';

export const _package = new DG.Package();

//name: PubChem
//tags: panel, widgets
//input: string molString {semType: Molecule}
//output: widget result
export async function pubChemPanel(molString: string): Promise<DG.Widget> {
  const pubChemId = await smilesToPubChem(molString);
  return new DG.Widget(ui.wait(async () => await init(pubChemId)));
}

//name: PubChem Substructure Search
//tags: panel, widgets
//input: string molString {semType: Molecule}
//output: widget result
export async function pubChemSubstructureSearchPanel(molString: string): Promise<DG.Widget> {
  return await pubChemSearchWidget(molString, 'substructure');
}

//name: PubChem Similarity Search
//tags: panel, widgets
//input: string molString {semType: Molecule}
//output: widget result
export async function pubChemSimilaritySearchPanel(molString: string): Promise<DG.Widget> {
  return await pubChemSearchWidget(molString, 'similarity');
}

//name: PubChem Identity Search
//tags: panel, widgets
//input: string molString {semType: Molecule}
//output: widget result
export async function pubChemIdentitySearch(molString: string): Promise<DG.Widget> {
  return await pubChemSearchWidget(molString, 'identity');
}

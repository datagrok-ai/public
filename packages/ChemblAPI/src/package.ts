/* eslint-disable no-unused-vars */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {chemblSearchWidget} from './utils';

export const _package = new DG.Package();
export * from './package.g';

export class PackageFunctions {
  @grok.decorators.panel({
    'name': 'Databases | ChEMBL | Substructure Search API',
    'meta': {role: 'widgets'},
  })
  static async chemblSubstructureSearchPanel(
    @grok.decorators.param({'options': {'semType': 'Molecule'}}) mol: string,
  ): Promise<DG.Widget> {
    return mol ? chemblSearchWidget(mol, true) : new DG.Widget(ui.divText('SMILES is empty'));
  }


  @grok.decorators.panel({
    'name': 'Databases | ChEMBL | Similarity Search API',
    'meta': {role: 'widgets'},
  })
  static async chemblSimilaritySearchPanel(
    @grok.decorators.param({'options': {'semType': 'Molecule'}}) mol: string,
  ): Promise<DG.Widget> {
    return mol ? chemblSearchWidget(mol) : new DG.Widget(ui.divText('SMILES is empty'));
  }


  @grok.decorators.func({'name': 'GetCompoundsIds', 'outputs': [{'type': 'object', 'name': 'result'}]})
  static async getCompoundsIds(inchiKey: string): Promise<{[key: string]: string | number}[]> {
    const url = `https://www.ebi.ac.uk/unichem/rest/inchikey/${inchiKey}`;
    const params: RequestInit = {method: 'GET', referrerPolicy: 'strict-origin-when-cross-origin'};
    const response = await grok.dapi.fetchProxy(url, params);
    const json = await response.json();
    return response.status !== 200 || json.error ? [] : json;
  }

  @grok.decorators.func({'name': 'Chembl Get by Id'})
  static async getById(id: string): Promise<DG.DataFrame> {
    if (!id.toLowerCase().startsWith('chembl'))
      id = `CHEMBL${id}`;

    try {
      return await grok.data.query(`${_package.name}:MoleculeJson`, {'molecule_chembl_id__exact': id});
    } catch (e: any) {
      grok.log.error(e.toString());
      return DG.DataFrame.create();
    }
  }
}

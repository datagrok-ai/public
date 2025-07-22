import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function chemblSearchWidget(mol: string, searchType: string): Promise<any> {
    return await grok.functions.call('ChEMBLAPI:ChemblSearchWidget', { mol, searchType });
  }

  export async function getById(): Promise<any> {
    return await grok.functions.call('ChEMBLAPI:GetById', {});
  }

  export async function chemblSubstructureSearchPanel(mol: string): Promise<any> {
    return await grok.functions.call('ChEMBLAPI:ChemblSubstructureSearchPanel', { mol });
  }

  export async function chemblSimilaritySearchPanel(mol: string): Promise<any> {
    return await grok.functions.call('ChEMBLAPI:ChemblSimilaritySearchPanel', { mol });
  }

  export async function getCompoundsIds(inchiKey: string): Promise<any> {
    return await grok.functions.call('ChEMBLAPI:GetCompoundsIds', { inchiKey });
  }
}

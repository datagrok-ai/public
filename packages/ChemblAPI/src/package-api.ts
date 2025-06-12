import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Funcs {
  export async function chemblSearchWidget(mol: string, searchType: string): Promise<any> {
    return await grok.functions.call('ChemblApi:ChemblSearchWidget', { mol, searchType });
  }

  export async function getById(): Promise<any> {
    return await grok.functions.call('ChemblApi:GetById', {});
  }

  export async function chemblSubstructureSearchPanel(mol: string): Promise<any> {
    return await grok.functions.call('ChemblApi:ChemblSubstructureSearchPanel', { mol });
  }

  export async function chemblSimilaritySearchPanel(mol: string): Promise<any> {
    return await grok.functions.call('ChemblApi:ChemblSimilaritySearchPanel', { mol });
  }

  export async function getCompoundsIds(inchiKey: string): Promise<any> {
    return await grok.functions.call('ChemblApi:GetCompoundsIds', { inchiKey });
  }
}

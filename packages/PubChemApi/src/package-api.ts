import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Funcs {
  export async function pubChemSubstructureSearchPanel(molString: string): Promise<any> {
    return await grok.functions.call('PubchemApi:PubChemSubstructureSearchPanel', { molString });
  }

  export async function pubChemSimilaritySearchPanel(molString: string): Promise<any> {
    return await grok.functions.call('PubchemApi:PubChemSimilaritySearchPanel', { molString });
  }

  export async function pubChemIdentitySearch(molString: string): Promise<any> {
    return await grok.functions.call('PubchemApi:PubChemIdentitySearch', { molString });
  }

  export async function pubChemToSmiles(id: string): Promise<any> {
    return await grok.functions.call('PubchemApi:PubChemToSmiles', { id });
  }

  export async function inchiKeysToSmiles(id: string): Promise<any> {
    return await grok.functions.call('PubchemApi:InchiKeysToSmiles', { id });
  }

  export async function getIupacName(smiles: string): Promise<any> {
    return await grok.functions.call('PubchemApi:GetIupacName', { smiles });
  }
}

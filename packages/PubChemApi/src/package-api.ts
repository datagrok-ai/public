import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function pubChemSubstructureSearchPanel(molString: string): Promise<any> {
    return await grok.functions.call('PubChem:PubChemSubstructureSearchPanel', { molString });
  }

  export async function pubChemSimilaritySearchPanel(molString: string): Promise<any> {
    return await grok.functions.call('PubChem:PubChemSimilaritySearchPanel', { molString });
  }

  export async function pubChemIdentitySearch(molString: string): Promise<any> {
    return await grok.functions.call('PubChem:PubChemIdentitySearch', { molString });
  }

  export async function pubChemToSmiles(id: string): Promise<any> {
    return await grok.functions.call('PubChem:PubChemToSmiles', { id });
  }

  export async function inchiKeysToSmiles(id: string): Promise<any> {
    return await grok.functions.call('PubChem:InchiKeysToSmiles', { id });
  }

  export async function getIupacName(smiles: string): Promise<any> {
    return await grok.functions.call('PubChem:GetIupacName', { smiles });
  }
}

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Queries {
  export async function searchPatentBySubstructure(pattern: string, maxMols: number): Promise<DG.DataFrame> {
    return await grok.data.query('Surechembl:SearchPatentBySubstructure', { pattern, maxMols });
  }

  export async function searchPatentBySimilarity(pattern: string, threshold: number, maxMols: number): Promise<DG.DataFrame> {
    return await grok.data.query('Surechembl:SearchPatentBySimilarity', { pattern, threshold, maxMols });
  }
}

export namespace Funcs {
  export async function sureChemblSubstructureSearchWidget(molecule: string): Promise<any> {
    return await grok.functions.call('Surechembl:SureChemblSubstructureSearchWidget', { molecule });
  }

  export async function sureChemblSimilaritySearchWidget(molecule: string): Promise<any> {
    return await grok.functions.call('Surechembl:SureChemblSimilaritySearchWidget', { molecule });
  }

  export async function sureChemblSubstructureSearch(molecule: string, limit: number): Promise<any> {
    return await grok.functions.call('Surechembl:SureChemblSubstructureSearch', { molecule, limit });
  }

  export async function sureChemblSimilaritySearch(molecule: string, limit: number, similarityThreshold: number): Promise<any> {
    return await grok.functions.call('Surechembl:SureChemblSimilaritySearch', { molecule, limit, similarityThreshold });
  }
}

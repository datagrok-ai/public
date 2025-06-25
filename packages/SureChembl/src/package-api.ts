import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace queries {
  export async function searchPatentBySubstructure(pattern: string, maxMols: number): Promise<DG.DataFrame> {
    return await grok.data.query('SureChEMBL:SearchPatentBySubstructure', { pattern, maxMols });
  }

  export async function searchPatentBySimilarity(pattern: string, threshold: number, maxMols: number): Promise<DG.DataFrame> {
    return await grok.data.query('SureChEMBL:SearchPatentBySimilarity', { pattern, threshold, maxMols });
  }
}

export namespace funcs {
  export async function sureChemblSubstructureSearchWidget(molecule: string): Promise<any> {
    return await grok.functions.call('SureChEMBL:SureChemblSubstructureSearchWidget', { molecule });
  }

  export async function sureChemblSimilaritySearchWidget(molecule: string): Promise<any> {
    return await grok.functions.call('SureChEMBL:SureChemblSimilaritySearchWidget', { molecule });
  }

  export async function sureChemblSubstructureSearch(molecule: string, limit: number): Promise<any> {
    return await grok.functions.call('SureChEMBL:SureChemblSubstructureSearch', { molecule, limit });
  }

  export async function sureChemblSimilaritySearch(molecule: string, limit: number, similarityThreshold: number): Promise<any> {
    return await grok.functions.call('SureChEMBL:SureChemblSimilaritySearch', { molecule, limit, similarityThreshold });
  }
}

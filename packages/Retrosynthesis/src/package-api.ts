import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function retroSynthesisPath(smiles: string): Promise<any> {
    return await grok.functions.call(':RetroSynthesisPath', { smiles });
  }

  export async function retrosynthesisTopMenu(): Promise<any> {
    return await grok.functions.call(':RetrosynthesisTopMenu', {});
  }

  //Generate retrosynthesis paths
  export async function retrosynthesisDemo(): Promise<any> {
    return await grok.functions.call(':RetrosynthesisDemo', {});
  }
}

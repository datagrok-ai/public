import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('Reinvent4:Info', {});
  }

  export async function getFolders(): Promise<any> {
    return await grok.functions.call('Reinvent4:GetFolders', {});
  }

  export async function reinventEditor(call: any): Promise<any> {
    return await grok.functions.call('Reinvent4:ReinventEditor', { call });
  }

  export async function runReinvent(ligand: string, optimize: string): Promise<any> {
    return await grok.functions.call('Reinvent4:RunReinvent', { ligand, optimize });
  }

  export async function reinvent(ligand: string, optimize: string): Promise<any> {
    return await grok.functions.call('Reinvent4:Reinvent', { ligand, optimize });
  }

  export async function reinventTopMenu(ligand: string, optimize: string): Promise<any> {
    return await grok.functions.call('Reinvent4:ReinventTopMenu', { ligand, optimize });
  }
}

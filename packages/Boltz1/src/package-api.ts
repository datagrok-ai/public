import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('Boltz1:Info', {});
  }

  export async function getBoltzConfigFolders(): Promise<any> {
    return await grok.functions.call('Boltz1:GetBoltzConfigFolders', {});
  }

  export async function runBoltz(config: string, msa: string): Promise<any> {
    return await grok.functions.call('Boltz1:RunBoltz', { config, msa });
  }

  export async function folding(df: DG.DataFrame, sequences: DG.Column): Promise<any> {
    return await grok.functions.call('Boltz1:Folding', { df, sequences });
  }

  export async function docking(df: DG.DataFrame, molecules: DG.Column, config: string): Promise<any> {
    return await grok.functions.call('Boltz1:Docking', { df, molecules, config });
  }

  export async function boltzWidget(molecule: any): Promise<any> {
    return await grok.functions.call('Boltz1:BoltzWidget', { molecule });
  }

  export async function isApplicableBoltz(molecule: string): Promise<any> {
    return await grok.functions.call('Boltz1:IsApplicableBoltz', { molecule });
  }

  export async function boltz1App(): Promise<any> {
    return await grok.functions.call('Boltz1:Boltz1App', {});
  }
}

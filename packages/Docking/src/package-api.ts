import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('Docking:Info', {});
  }

  export async function getAutoDockService(): Promise<any> {
    return await grok.functions.call('Docking:GetAutoDockService', {});
  }

  export async function autoDockApp(): Promise<any> {
    return await grok.functions.call('Docking:AutoDockApp', {});
  }

  export async function getConfigFiles(): Promise<any> {
    return await grok.functions.call('Docking:GetConfigFiles', {});
  }

  export async function dockLigandCached(jsonForm: string, containerId: string): Promise<any> {
    return await grok.functions.call('Docking:DockLigandCached', { jsonForm, containerId });
  }

  //Autodock plugin UI
  export async function runAutodock5(table: DG.DataFrame, ligands: DG.Column, target: string, poses: number): Promise<any> {
    return await grok.functions.call('Docking:RunAutodock5', { table, ligands, target, poses });
  }

  export async function isApplicableAutodock(molecule: string): Promise<any> {
    return await grok.functions.call('Docking:IsApplicableAutodock', { molecule });
  }

  export async function autodockWidget(molecule: any): Promise<any> {
    return await grok.functions.call('Docking:AutodockWidget', { molecule });
  }

  export async function getAutodockSingle(): Promise<any> {
    return await grok.functions.call('Docking:GetAutodockSingle', {});
  }

  //Small molecule docking to a macromolecule with pose visualization
  export async function demoDocking(): Promise<any> {
    return await grok.functions.call('Docking:DemoDocking', {});
  }

  export async function autodockPanel(smiles: any): Promise<any> {
    return await grok.functions.call('Docking:AutodockPanel', { smiles });
  }

  export async function dockingView(path: string): Promise<any> {
    return await grok.functions.call('Docking:DockingView', { path });
  }
}

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('@datagrok/docking:Info', {});
  }

  export async function getAutoDockService(): Promise<any> {
    return await grok.functions.call('@datagrok/docking:GetAutoDockService', {});
  }

  export async function autoDockApp(): Promise<any> {
    return await grok.functions.call('@datagrok/docking:AutoDockApp', {});
  }

  export async function getConfigFiles(): Promise<any> {
    return await grok.functions.call('@datagrok/docking:GetConfigFiles', {});
  }

  export async function dockLigandCached(jsonForm: string, containerId: string): Promise<any> {
    return await grok.functions.call('@datagrok/docking:DockLigandCached', { jsonForm, containerId });
  }

  //Autodock plugin UI
  export async function runAutodock5(table: DG.DataFrame, ligands: DG.Column, target: string, poses: number): Promise<any> {
    return await grok.functions.call('@datagrok/docking:RunAutodock5', { table, ligands, target, poses });
  }

  export async function isApplicableAutodock(molecule: string): Promise<any> {
    return await grok.functions.call('@datagrok/docking:IsApplicableAutodock', { molecule });
  }

  export async function autodockWidget(molecule: any): Promise<any> {
    return await grok.functions.call('@datagrok/docking:AutodockWidget', { molecule });
  }

  export async function getAutodockSingle(): Promise<any> {
    return await grok.functions.call('@datagrok/docking:GetAutodockSingle', {});
  }

  //Small molecule docking to a macromolecule with pose visualization
  export async function demoDocking(): Promise<any> {
    return await grok.functions.call('@datagrok/docking:DemoDocking', {});
  }

  export async function autodockPanel(smiles: any): Promise<any> {
    return await grok.functions.call('@datagrok/docking:AutodockPanel', { smiles });
  }

  export async function dockingView(path: string): Promise<any> {
    return await grok.functions.call('@datagrok/docking:DockingView', { path });
  }
}

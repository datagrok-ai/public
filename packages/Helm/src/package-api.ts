import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function initHelm(): Promise<any> {
    return await grok.functions.call('@datagrok/helm:InitHelm', {});
  }

  //Helm renderer service
  export async function getHelmService(): Promise<any> {
    return await grok.functions.call('@datagrok/helm:GetHelmService', {});
  }

  export async function helmCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/helm:HelmCellRenderer', {});
  }

  //Macromolecule
  export async function editMoleculeCell(cell: any): Promise<any> {
    return await grok.functions.call('@datagrok/helm:EditMoleculeCell', { cell });
  }

  //Adds editor
  export async function openEditor(mol: any): Promise<any> {
    return await grok.functions.call('@datagrok/helm:OpenEditor', { mol });
  }

  export async function propertiesWidget(sequence: any): Promise<any> {
    return await grok.functions.call('@datagrok/helm:PropertiesWidget', { sequence });
  }

  export async function getMolfiles(col: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/helm:GetMolfiles', { col });
  }

  export async function helmInput(name: string, options: any): Promise<any> {
    return await grok.functions.call('@datagrok/helm:HelmInput', { name, options });
  }

  export async function getHelmHelper(): Promise<any> {
    return await grok.functions.call('@datagrok/helm:GetHelmHelper', {});
  }

  export async function measureCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/helm:MeasureCellRenderer', {});
  }

  export async function highlightMonomers(): Promise<any> {
    return await grok.functions.call('@datagrok/helm:HighlightMonomers', {});
  }
}

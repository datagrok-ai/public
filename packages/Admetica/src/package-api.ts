import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('@datagrok/admetica:Info', {});
  }

  export async function runAdmetica(csvString: string, queryParams: string, raiseException: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/admetica:RunAdmetica', { csvString, queryParams, raiseException });
  }

  export async function admeticaWidget(smiles: any): Promise<any> {
    return await grok.functions.call('@datagrok/admetica:AdmeticaWidget', { smiles });
  }

  export async function getModels(property: string): Promise<any> {
    return await grok.functions.call('@datagrok/admetica:GetModels', { property });
  }

  export async function admeticaHT(table: DG.DataFrame, molecules: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/admetica:AdmeticaHT', { table, molecules });
  }

  export async function admeticaEditor(call: any): Promise<any> {
    return await grok.functions.call('@datagrok/admetica:AdmeticaEditor', { call });
  }

  export async function admeticaMenu(table: DG.DataFrame, molecules: DG.Column, template: string, addPiechart: boolean, addForm: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/admetica:AdmeticaMenu', { table, molecules, template, addPiechart, addForm });
  }

  export async function admeProperty(molecule: string, prop: string): Promise<any> {
    return await grok.functions.call('@datagrok/admetica:AdmeProperty', { molecule, prop });
  }

  export async function admeticaApp(): Promise<any> {
    return await grok.functions.call('@datagrok/admetica:AdmeticaApp', {});
  }

  //Evaluating ADMET properties
  export async function admeticaDemo(): Promise<any> {
    return await grok.functions.call('@datagrok/admetica:AdmeticaDemo', {});
  }
}

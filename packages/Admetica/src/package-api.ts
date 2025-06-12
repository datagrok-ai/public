import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('Admetica:Info', {});
  }

  export async function runAdmetica(csvString: string, queryParams: string, raiseException: boolean): Promise<any> {
    return await grok.functions.call('Admetica:RunAdmetica', { csvString, queryParams, raiseException });
  }

  export async function admeticaWidget(smiles: any): Promise<any> {
    return await grok.functions.call('Admetica:AdmeticaWidget', { smiles });
  }

  export async function getModels(property: string): Promise<any> {
    return await grok.functions.call('Admetica:GetModels', { property });
  }

  export async function admeticaHT(table: DG.DataFrame, molecules: DG.Column): Promise<any> {
    return await grok.functions.call('Admetica:AdmeticaHT', { table, molecules });
  }

  export async function admeticaEditor(call: any): Promise<any> {
    return await grok.functions.call('Admetica:AdmeticaEditor', { call });
  }

  export async function admeticaMenu(table: DG.DataFrame, molecules: DG.Column, template: string, addPiechart: boolean, addForm: boolean): Promise<any> {
    return await grok.functions.call('Admetica:AdmeticaMenu', { table, molecules, template, addPiechart, addForm });
  }

  export async function admeProperty(molecule: string, prop: string): Promise<any> {
    return await grok.functions.call('Admetica:AdmeProperty', { molecule, prop });
  }

  export async function admeticaApp(): Promise<any> {
    return await grok.functions.call('Admetica:AdmeticaApp', {});
  }

  //Evaluating ADMET properties
  export async function admeticaDemo(): Promise<any> {
    return await grok.functions.call('Admetica:AdmeticaDemo', {});
  }
}

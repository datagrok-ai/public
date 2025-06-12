import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Funcs {
  export async function jdxFileHandler(bytes: string): Promise<any> {
    return await grok.functions.call('Nmrium:JdxFileHandler', { bytes });
  }

  export async function dxFileHandler(bytes: string): Promise<any> {
    return await grok.functions.call('Nmrium:DxFileHandler', { bytes });
  }

  export async function nmriumFileHandler(bytes: string): Promise<any> {
    return await grok.functions.call('Nmrium:NmriumFileHandler', { bytes });
  }

  export async function previewNMRData(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('Nmrium:PreviewNMRData', { file });
  }

  export async function previewNMRFromDX(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('Nmrium:PreviewNMRFromDX', { file });
  }

  export async function checkNmriumJdx(content: string): Promise<any> {
    return await grok.functions.call('Nmrium:CheckNmriumJdx', { content });
  }
}

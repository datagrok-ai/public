import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function jdxFileHandler(bytes: string): Promise<any> {
    return await grok.functions.call('@datagrok/nmrium:JdxFileHandler', { bytes });
  }

  export async function dxFileHandler(bytes: string): Promise<any> {
    return await grok.functions.call('@datagrok/nmrium:DxFileHandler', { bytes });
  }

  export async function nmriumFileHandler(bytes: string): Promise<any> {
    return await grok.functions.call('@datagrok/nmrium:NmriumFileHandler', { bytes });
  }

  export async function previewNMRData(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/nmrium:PreviewNMRData', { file });
  }

  export async function previewNMRFromDX(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/nmrium:PreviewNMRFromDX', { file });
  }

  export async function checkNmriumJdx(content: string): Promise<any> {
    return await grok.functions.call('@datagrok/nmrium:CheckNmriumJdx', { content });
  }
}

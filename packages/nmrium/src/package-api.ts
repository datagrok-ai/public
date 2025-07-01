import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function jdxFileHandler(bytes: string ): Promise<any> {
    return await grok.functions.call('NMRium:JdxFileHandler', { bytes });
  }

  export async function dxFileHandler(bytes: string ): Promise<any> {
    return await grok.functions.call('NMRium:DxFileHandler', { bytes });
  }

  export async function nmriumFileHandler(bytes: string ): Promise<any> {
    return await grok.functions.call('NMRium:NmriumFileHandler', { bytes });
  }

  export async function previewNMRData(file: DG.FileInfo ): Promise<DG.View> {
    return await grok.functions.call('NMRium:PreviewNMRData', { file });
  }

  export async function previewNMRFromDX(file: DG.FileInfo ): Promise<DG.View> {
    return await grok.functions.call('NMRium:PreviewNMRFromDX', { file });
  }

  export async function checkNmriumJdx(content: string ): Promise<boolean> {
    return await grok.functions.call('NMRium:CheckNmriumJdx', { content });
  }
}

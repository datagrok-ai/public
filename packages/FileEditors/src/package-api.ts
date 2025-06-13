import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('@datagrok/file-editors:Info', {});
  }

  export async function previewPdf(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/file-editors:PreviewPdf', { file });
  }

  export async function viewPdf(bytes: any): Promise<any> {
    return await grok.functions.call('@datagrok/file-editors:ViewPdf', { bytes });
  }

  export async function previewDocx(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/file-editors:PreviewDocx', { file });
  }

  export async function viewDocx(bytes: any): Promise<any> {
    return await grok.functions.call('@datagrok/file-editors:ViewDocx', { bytes });
  }

  export async function previewRtf(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/file-editors:PreviewRtf', { file });
  }
}

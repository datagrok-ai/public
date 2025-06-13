import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('@datagrok/genome-browser:Info', {});
  }

  export async function previewGenomeFileBrowse(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/genome-browser:PreviewGenomeFileBrowse', { file });
  }
}

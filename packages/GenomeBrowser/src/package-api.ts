import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('GenomeBrowser:Info', {});
  }

  export async function previewGenomeFileBrowse(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('GenomeBrowser:PreviewGenomeFileBrowse', { file });
  }
}

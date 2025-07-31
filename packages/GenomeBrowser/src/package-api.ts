import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function info(): Promise<void> {
    return await grok.functions.call('GenomeBrowser:Info', {});
  }

  export async function previewGenomeFileBrowse(file: DG.FileInfo ): Promise<DG.View> {
    return await grok.functions.call('GenomeBrowser:PreviewGenomeFileBrowse', { file });
  }
}

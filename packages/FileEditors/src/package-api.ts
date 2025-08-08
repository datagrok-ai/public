import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function info(): Promise<void> {
    return await grok.functions.call('FileEditors:Info', {});
  }

  export async function previewPdf(file: DG.FileInfo ): Promise<DG.View> {
    return await grok.functions.call('FileEditors:PreviewPdf', { file });
  }

  export async function viewPdf(bytes: any ): Promise<any> {
    return await grok.functions.call('FileEditors:ViewPdf', { bytes });
  }

  export async function previewDocx(file: DG.FileInfo ): Promise<DG.View> {
    return await grok.functions.call('FileEditors:PreviewDocx', { file });
  }

  export async function viewDocx(bytes: any ): Promise<any> {
    return await grok.functions.call('FileEditors:ViewDocx', { bytes });
  }

  export async function previewRtf(file: DG.FileInfo ): Promise<DG.View> {
    return await grok.functions.call('FileEditors:PreviewRtf', { file });
  }

  export async function previewTex(file: DG.FileInfo ): Promise<DG.View> {
    return await grok.functions.call('FileEditors:PreviewTex', { file });
  }
}

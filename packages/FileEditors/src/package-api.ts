import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('FileEditors:Info', {});
  }

  export async function previewPdf(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('FileEditors:PreviewPdf', { file });
  }

  export async function viewPdf(bytes: any): Promise<any> {
    return await grok.functions.call('FileEditors:ViewPdf', { bytes });
  }

  export async function previewDocx(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('FileEditors:PreviewDocx', { file });
  }

  export async function viewDocx(bytes: any): Promise<any> {
    return await grok.functions.call('FileEditors:ViewDocx', { bytes });
  }

  export async function previewRtf(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('FileEditors:PreviewRtf', { file });
  }
}

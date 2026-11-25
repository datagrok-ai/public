/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../css/file_editors.css';
import {getDocument, PDFDocumentProxy, PDFPageProxy} from 'pdfjs-dist';
import 'pdfjs-dist/webpack';
import {LatexViewer} from './latex/latex-viewer';
import {renderOdt} from './odt';
import {readDocParagraphs} from './doc';
import {isSafeUrl, renderDocx} from './docx';
export * from './package.g';
export const _package = new DG.Package();

/** A view whose content scrolls as one document; a render error is shown in place of the content. */
async function documentView(name: string, render: (host: HTMLDivElement) => Promise<unknown>): Promise<DG.View> {
  const view = DG.View.create();
  view.name = name;
  const host = ui.div([], 'fe-document-preview');
  view.append(host);
  try {
    await render(host);
  }
  catch (e) {
    host.replaceChildren(ui.divText(`${name} cannot be displayed: ${e instanceof Error ? e.message : e}`));
  }
  return view;
}


export class PackageFunctions {
  @grok.decorators.func()
  static info() {
    grok.shell.info(_package.webRoot);
  }


  @grok.decorators.fileViewer({fileViewer: 'pdf'})
  static async previewPdf(
    file: DG.FileInfo): Promise<DG.View> {
    const view = DG.View.create();
    view.name = file.name;
    const root = ui.div();
    root.classList.add('fe-pdf-preview');
    view.append(root);
    let currPage = 1;
    let numPages = 0;
    let pdf: PDFDocumentProxy | null = null;
    const bytes = await file.readAsBytes();
    getDocument(bytes).promise.then((pdf_) => {
      pdf = pdf_;
      numPages = pdf_.numPages >= 10 ? 10 : pdf_.numPages;
      pdf.getPage(1).then(handlePage);
    });
    return view;

    function handlePage(page: PDFPageProxy) {
      const canvas = ui.canvas(470);
      root.appendChild(canvas);
      const scale = 0.6;
      const viewport = page.getViewport({scale: canvas.width / page.getViewport({scale: scale}).width});
      const context = canvas.getContext('2d')!;
      canvas.width = 470;
      canvas.height = viewport.height;
      page.render({
        canvasContext: context,
        viewport: viewport,
      });
      currPage++;
      if (pdf && currPage <= numPages)
        pdf.getPage(currPage).then(handlePage);
    }
  }


  @grok.decorators.fileHandler({ext: 'pdf'})
  static async viewPdf(
    @grok.decorators.param({type: 'list'}) bytes: Uint8Array): Promise<DG.DataFrame[]> {
    bytes = new Uint8Array(bytes);
    const view = DG.View.create();
    view.name = 'PDF File';
    const root = ui.div();
    root.classList.add('pv-pdf-preview');
    view.append(root);
    let currPage = 1;
    let pdf: PDFDocumentProxy | null = null;
    getDocument(bytes).promise.then((pdf_) => {
      pdf = pdf_;
      pdf.getPage(1).then(handlePage);
    });
    grok.shell.addView(view);
    return [];

    function handlePage(page: PDFPageProxy) {
      const canvas = ui.canvas();
      root.appendChild(canvas);
      const scale = 0.25;
      const viewport = page.getViewport({scale: canvas.width / page.getViewport({scale: scale}).width});
      const context = canvas.getContext('2d')!;
      canvas.width = viewport.width;
      canvas.height = viewport.height;
      page.render({
        canvasContext: context,
        viewport: viewport,
      });
      currPage++;
      if (pdf && currPage <= pdf.numPages)
        pdf.getPage(currPage).then(handlePage);
    }
  }


  @grok.decorators.fileViewer({fileViewer: 'docx'})
  static async previewDocx(
    file: DG.FileInfo): Promise<DG.View> {
    return documentView(file.name, async (host) => renderDocx(await file.readAsBytes(), host));
  }


  @grok.decorators.fileHandler({ext: 'docx'})
  static async viewDocx(
    @grok.decorators.param({type: 'list'}) bytes: Uint8Array): Promise<DG.DataFrame[]> {
    grok.shell.addView(await documentView('DOCX File', (host) => renderDocx(new Uint8Array(bytes), host)));
    return [];
  }


  @grok.decorators.fileViewer({fileViewer: 'rtf'})
  static async previewRtf(
    file: DG.FileInfo): Promise<DG.View> {
    return documentView(file.name, async (host) => {
      const {RTFJS, EMFJS, WMFJS} = await import('rtf.js');
      for (const lib of [RTFJS, EMFJS, WMFJS])
        lib.loggingEnabled(false);
      const elements = await new RTFJS.Document(new Uint8Array(await file.readAsBytes()).buffer, {
        onHyperlink: (create, link) => {
          const a = create();
          if (!isSafeUrl(link.url()))
            a.removeAttribute('href');
          return {element: a, content: a};
        },
      }).render();
      for (const element of elements) {
        element.style.minHeight = '';
        host.append(element);
      }
    });
  }


  @grok.decorators.fileViewer({fileViewer: 'odt'})
  static async previewOdt(
    file: DG.FileInfo): Promise<DG.View> {
    return documentView(file.name, async (host) => renderOdt(await file.readAsBytes(), host));
  }


  @grok.decorators.fileViewer({fileViewer: 'doc'})
  static async previewDoc(
    file: DG.FileInfo): Promise<DG.View> {
    return documentView(file.name, async (host) => {
      for (const paragraph of readDocParagraphs(await file.readAsBytes()))
        host.append(ui.p(paragraph || '\u00a0'));
    });
  }


  @grok.decorators.fileViewer({fileViewer: 'tex'})
  static async previewTex(
    file: DG.FileInfo): Promise<DG.View> {
    // const view = DG.View.create();
    // view.name = file.fileName;

    // try {
    //   const latexText = await file.readAsString();
    //   view.append(getElementWithLatexContent(latexText));
    // } catch (err) {
    //   view.append(ui.h2('The file is corrupted and cannot be opened!'));
    //   if (err instanceof Error)
    //     grok.shell.error(err.message);
    // }

    // return view;

    const viewer = await LatexViewer.create(file);

    return viewer.view;
  }
}

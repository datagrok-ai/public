/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../css/file_editors.css';
import {getDocument, PDFDocumentProxy, PDFPageProxy} from 'pdfjs-dist';
import 'pdfjs-dist/webpack';
import {renderAsync} from 'docx-preview';
import {RTFJS} from "rtf.js";

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: fileViewer
//meta.fileViewer: pdf
//input: file file
//output: view v
export async function previewPdf(file: DG.FileInfo): Promise<DG.View> {
  const view = DG.View.create();
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

//tags: file-handler
//meta.ext: pdf
//input: list bytes
//output: list tables
export async function viewPdf(bytes: Uint8Array): Promise<DG.DataFrame[]> {
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

//tags: fileViewer
//meta.fileViewer: docx
//input: file file
//output: view v
export async function previewDocx(file: DG.FileInfo): Promise<DG.View> {
  const view = DG.View.create();
  view.root.classList.add('fe-docx-preview');
  const bytes = await file.readAsBytes();
  renderAsync(bytes, view.root, undefined, {inWrapper: false, ignoreWidth: true, ignoreHeight: true});
  return view;
}

//tags: file-handler
//meta.ext: docx
//input: list bytes
//output: list tables
export async function viewDocx(bytes: Uint8Array): Promise<DG.DataFrame[]> {
  const view = DG.View.create();
  bytes = new Uint8Array(bytes);
  view.name = 'DOCX File';
  renderAsync(bytes, view.root, undefined, {inWrapper: false});
  grok.shell.addView(view);
  return [];
}

//tags: fileViewer
//meta.fileViewer: rtf
//input: file file
//output: view v
export async function previewRtf(file: DG.FileInfo): Promise<DG.View> {
  const view = DG.View.create();

  file.readAsBytes().then((bytes) => {
    const doc = new RTFJS.Document(bytes, {});
    doc.render().then((elements) => {
      view.root.append(...elements);
    })

  });

  return view;
}

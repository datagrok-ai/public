import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// @ts-ignore
import {parse, HtmlGenerator} from 'latex.js';
import {_package} from './package';

export function getElementWithLatexContent(latexText: string): HTMLIFrameElement {
  const generator = new HtmlGenerator({hyphenate: false});
  const htmlGenerator = parse(latexText, {generator: generator});
  const doc = htmlGenerator.htmlDocument(`${_package.webRoot}css/latex/`) as Document;

  const iframe = ui.iframe({width: '100%', height: '100%'});
  const rightMargin = doc.body.children[1];

  if ((rightMargin !== null) && (rightMargin !== null))
    (rightMargin as HTMLElement).style.display = 'none';

  iframe.srcdoc = doc.documentElement.outerHTML;
  iframe.style.border = 'transparent';

  return iframe;
}

/** */
export class LatexViewer {
  public view: DG.View;
  private file: DG.FileInfo;
  private isPlatformFile: boolean = false;

  static async create(file: DG.FileInfo): Promise<LatexViewer> {
    try {
      const latexText = await file.readAsString();
      const exist = await grok.dapi.files.exists(file);

      console.log('Platform file: ', exist);

      return new LatexViewer(file, latexText, exist);
    } catch (err) {
      if (err instanceof Error)
        grok.shell.error(err.message);

      return new LatexViewer(file, null, false);
    }
  }

  private constructor(file: DG.FileInfo, latexText: string | null, exist: boolean) {
    this.file = file;
    this.view = DG.View.create();
    this.view.name = file.fileName;
    this.isPlatformFile = exist;

    if (latexText === null)
      this.view.append(ui.h2('The file is corrupted and cannot be opened!'));
    else {
      try {
        this.view.append(getElementWithLatexContent(latexText));
      } catch (err) {
        if (err instanceof Error)
          grok.shell.error(err.message);

        this.view.append(ui.h2('LaTeX code contains errors!'));
      }

      this.buildIO();
    }
  };

  private buildIO(): void {}
}; // LatexViewer

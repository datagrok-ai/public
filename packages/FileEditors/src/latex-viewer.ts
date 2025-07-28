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

  private codeDiv = ui.div(['Code']);
  private contentDiv = ui.div();
  private container = ui.div();

  private isEditorShown = false;

  static async create(file: DG.FileInfo): Promise<LatexViewer> {
    try {
      const latexText = await file.readAsString();
      const exist = await grok.dapi.files.exists(file);

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
        this.contentDiv.append(getElementWithLatexContent(latexText));
      } catch (err) {
        if (err instanceof Error)
          grok.shell.error(err.message);

        this.view.append(ui.h2('LaTeX code contains errors!'));
      }

      this.buildIO(latexText);
    }
  };

  private buildIO(latexText: string): void {
    this.container = ui.divH([this.codeDiv, this.contentDiv], {style: {width: '100%'}});
    this.updateEditorVisibility(this.isEditorShown);

    this.view.append(this.container);
    this.view.setRibbonPanels([[
      this.getEditToggle(),
    ]]);
  }

  private getEditToggle(): HTMLElement {
    const input = ui.input.toggle('', {value: this.isEditorShown});

    const span = ui.span(['Edit']);
    span.classList.add('fe-latex-viewer-ribbon-text');

    const wgt = ui.divH([input.root, span]);

    ui.tooltip.bind(wgt, this.isEditorShown ? 'Finish editing' : 'Edit source');

    wgt.onclick = (e) => {
      e.stopImmediatePropagation();
      e.preventDefault();

      this.isEditorShown = !this.isEditorShown;
      input.value = this.isEditorShown;
      this.updateEditorVisibility(this.isEditorShown);
      ui.tooltip.bind(wgt, this.isEditorShown ? 'Finish editing' : 'Edit source');
      // this.tabControl.currentPane = this.isEditState ? this.editPane : this.solvePane;
      // this.updateRibbonWgts();
      // this.updateRefreshWidget(this.isModelChanged);
    };

    return wgt;
  } // getEditToggle

  private updateEditorVisibility(toShow: boolean): void {
    this.codeDiv.hidden = !toShow;
    this.contentDiv.style.width = toShow ? '50%' : '100%';
    this.codeDiv.style.width = toShow ? '50%' : '0%';
  }
}; // LatexViewer

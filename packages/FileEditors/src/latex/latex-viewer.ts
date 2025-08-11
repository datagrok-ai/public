/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// @ts-ignore
import {parse, HtmlGenerator} from 'latex.js';
import {_package} from '../package';

import {basicSetup, EditorView} from 'codemirror';
import {latex} from 'codemirror-lang-latex';
import {HIGH, LIMIT, LOW, MARGIN_IDX} from './constants';
import {debounce, isActionKey, isCutPaste} from './utils';

/** Viewer for tex-files */
export class LatexViewer {
  public view: DG.View;
  private file: DG.FileInfo;
  private isPlatformFile = false;
  private debounceMs = HIGH;

  private codeDiv = ui.div();
  private contentDiv = ui.div();
  private container = ui.div();

  private isEditorShown = false;
  private isLatexCodeChanged = false;

  private editorView: EditorView;

  private saveIcn = ui.iconFA('save', async () => {
    if (this.isLatexCodeChanged) {
      await this.saveFile();
      this.isLatexCodeChanged = false;
    }
  });

  private prevNode: Node | null = null;

  /** Load tex-file and create the viewer */
  static async create(file: DG.FileInfo): Promise<LatexViewer> {
    try {
      const latexText = await file.readAsString();
      console.log('Length:', latexText.length);
      const exist = await grok.dapi.files.exists(file);

      return new LatexViewer(file, latexText, exist);
    } catch (err) {
      if (err instanceof Error)
        grok.shell.error(err.message);

      return new LatexViewer(file, null, false);
    }
  }

  private updateSaveWgt(enabled: boolean): void {
    this.saveIcn.style.color = this.getColor(enabled);
    ui.tooltip.bind(
      this.saveIcn,
      enabled ? (this.isPlatformFile ? 'Save changes' : 'Save to My files') : 'No unsaved changes',
    );
  }

  private updateDebounceMs(text: string): void {
    this.debounceMs = text.length > LIMIT ? HIGH : LOW;
  }

  private constructor(file: DG.FileInfo, latexText: string | null, exist: boolean) {
    this.file = file;
    this.view = DG.View.create();
    this.view.name = file.fileName;
    this.isPlatformFile = exist;

    if (latexText === null)
      this.view.append(ui.h2('The file is corrupted and cannot be opened!'));
    else {
      this.updateDebounceMs(latexText);

      try {
        this.prevNode = this.contentDiv.appendChild(this.getElementWithLatexContent(latexText));
      } catch (err) {
        if (err instanceof Error)
          grok.shell.error(err.message);

        this.prevNode = this.contentDiv.appendChild(ui.h2('LaTeX code contains errors!'));
      }

      this.editorView = new EditorView({
        doc: latexText,
        extensions: [basicSetup, latex()],
      });

      this.buildIO();
    }
  };

  private getElementWithLatexContent(latexText: string): HTMLIFrameElement {
    const generator = new HtmlGenerator({hyphenate: false});
    const htmlGenerator = parse(latexText, {generator: generator});
    const doc = htmlGenerator.htmlDocument(`${_package.webRoot}css/latex/`) as Document;

    const iframe = ui.iframe({width: '100%', height: '100%'});
    const rightMargin = doc.body.children[MARGIN_IDX];

    if ((rightMargin !== null) && (rightMargin !== null))
      (rightMargin as HTMLElement).style.display = 'none';

    iframe.srcdoc = doc.documentElement.outerHTML;
    iframe.style.border = 'transparent';

    return iframe;
  }

  private buildIO(): void {
    this.codeDiv.append(this.editorView.dom);
    this.codeDiv.style.height = '100%';
    this.contentDiv.style.height = '100%';
    this.container = ui.divH([this.codeDiv, this.contentDiv], {style: {width: '100%', height: '100%'}});
    this.updateEditorVisibility(this.isEditorShown);

    const handleKeyPress = (event: KeyboardEvent) => {
      if (isActionKey(event.key) || isCutPaste(event)) {
        this.isLatexCodeChanged = true;
        this.updateSaveWgt(true);
        this.applyChanges(this.editorView.state.doc.toString());
      }
    };

    const debouncedInput = debounce(handleKeyPress, this.debounceMs);
    this.editorView.dom.addEventListener('keydown', debouncedInput);
    this.updateSaveWgt(false);

    this.view.append(this.container);
    this.view.setRibbonPanels([[
      this.saveIcn,
      this.getDownLoadIcon(),
      this.getEditToggle(),
    ]]);
  }

  private async saveFile(): Promise<void> {
    if (this.isPlatformFile) {
      try {
        const source = new DG.FileSource();

        await source.writeAsText(this.file, this.editorView!.state.doc.toString());
        grok.shell.info('Saved');
        this.updateSaveWgt(false);
      } catch (err) {
        grok.shell.error(`Failed to save changes: ${err instanceof Error ? err.message : 'the platform issue'}`);
      }

      return;
    }

    const texCode = this.editorView!.state.doc.toString();
    let fileName = this.file.name.replaceAll(' ', '-');
    const login = grok.shell.user.project.name;
    const folder = login.charAt(0).toUpperCase() + login.slice(1) + ':Home/';
    const files = await grok.dapi.files.list(folder);
    const existingNames = files.filter((file) => file.extension === 'tex').map((file) => file.name);

    const nameInput = ui.input.string('Name', {
      value: fileName,
      nullable: false,
      onValueChanged: () => {
        if (nameInput.value !== null) {
          fileName = nameInput.value.replaceAll(' ', '-');
          dlg.getButton('Save').disabled = false;
        } else
          dlg.getButton('Save').disabled = true;
      },
    });

    const save = async () => {
      const path = `${folder}${fileName}`;

      try {
        await grok.dapi.files.writeAsText(path, texCode);
        grok.shell.info(`Saved to ${path}`);
        this.updateSaveWgt(false);
      } catch (error) {
        grok.shell.error(`Failed to save: ${error instanceof Error ? error.message : 'platform issue'}`);
      }

      dlg.close();
    };

    const dlg = ui.dialog({title: 'Save to My files'})
      .add(nameInput)
      .addButton('Save', async () => {
        if (!existingNames.includes(fileName))
          await save();
        else {
          ui.dialog({title: 'WARNING'})
            .add(ui.label('Overwrite existing file?'))
            .onOK(async () => await save())
            .show();
        }
      }, undefined, 'Save tex-file to My files')
      .show();
  }; // saveToMyFiles

  private getDownLoadIcon(): HTMLElement {
    const icon = ui.iconFA('arrow-to-bottom', () => {
      const link = document.createElement('a');
      const file = new Blob([this.editorView!.state.doc.toString()], {type: 'text/plain'});
      link.href = URL.createObjectURL(file);
      link.download = this.file.name;
      link.click();
      URL.revokeObjectURL(link.href);
    }, 'Download tex-file');
    icon.classList.add('fe-latex-viewer-ribbon-download');

    return icon;
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
    };

    return wgt;
  } // getEditToggle

  private updateEditorVisibility(toShow: boolean): void {
    this.codeDiv.hidden = !toShow;
    this.contentDiv.style.width = toShow ? '50%' : '100%';
    this.codeDiv.style.width = toShow ? '50%' : '0%';
    this.contentDiv.style.borderLeft = toShow ? '1px solid var(--steel-1)' : 'none';
  }

  private getColor(enabled: boolean) {
    return enabled ? '#40607F' : 'var(--grey-3)';
  }

  private applyChanges(latexText: string): void {
    const scrollTop = this.getScrollTop();

    if (this.prevNode !== null)
      this.contentDiv.removeChild(this.prevNode);

    try {
      const newIframe = this.getElementWithLatexContent(latexText);
      this.prevNode = this.contentDiv.appendChild(newIframe);
      newIframe.style.opacity = '0';
      newIframe.style.transition = 'opacity 300ms ease-in';
      newIframe.style.opacity = '1';

      this.updateDebounceMs(latexText);

      newIframe.addEventListener('load', () => {
        requestAnimationFrame(() => {
          const maxScrollTop = newIframe.contentDocument.documentElement.scrollHeight -
          newIframe.contentDocument.documentElement.clientHeight;
          newIframe.contentDocument.documentElement.scrollTop = Math.min(scrollTop, Math.max(0, maxScrollTop));
        });
      });
    } catch (err) {
      if (err instanceof Error)
        grok.shell.error(err.message);

      this.prevNode = this.contentDiv.appendChild(ui.h2('LaTeX code contains errors!'));
    }
  }

  private getScrollTop(): number {
    const iframe = this.contentDiv.querySelector('iframe');

    if ((iframe === null) || (iframe === undefined))
      return 0;

    const scrollTop = iframe.contentDocument?.documentElement?.scrollTop ?? 0;
    iframe.style.transition = 'opacity 300ms ease-out';
    iframe.style.opacity = '0';

    return scrollTop;
  } // applyChanges
}; // LatexViewer

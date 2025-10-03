import * as ui from 'datagrok-api/ui';
import {UiUtils} from '@datagrok-libraries/compute-utils';

export class FileInputUtils {
  static createFileInputPane(fileCallback: (file: File) => Promise<void>): HTMLElement {
    const fileInputEditor = this.initializeFileInputEditor(fileCallback);
    this.removeLabels(fileInputEditor.root);
    this.styleInputEditor(fileInputEditor.root);
    this.setupDragAndDrop(fileInputEditor.root);
    this.removeOptionsIcon(fileInputEditor.root);

    fileInputEditor.root.classList.add('demo-file-input');
    return ui.divV([fileInputEditor], {classes: 'demo-file-input-container'});
  }

  private static initializeFileInputEditor(fileCallback: (file: File) => Promise<void>) {
    const fileInputEditor = UiUtils.fileInput('', null, async (file: File) => {
      await fileCallback(file);
    }, null);

    fileInputEditor.stringValue = 'Drag and drop a CSV file here, or click to select a file.';
    return fileInputEditor;
  }

  private static removeLabels(root: HTMLElement): void {
    const labels = root.querySelectorAll<HTMLElement>('.ui-label, .ui-input-label');
    labels.forEach((label) => label.remove());
  }

  private static styleInputEditor(root: HTMLElement): void {
    const inputEditor = root.querySelector<HTMLElement>('.ui-input-editor');
    if (!inputEditor) return;

    Object.assign(inputEditor.style, {
      width: '100%',
      height: '100%',
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      backgroundColor: '#ffffff',
      color: '#007bff',
      fontSize: '14px',
      cursor: 'pointer',
      textAlign: 'center',
      borderBottom: 'none',
      border: '1px dashed #007bff',
    });
  }

  private static setupDragAndDrop(root: HTMLElement): void {
    const inputEditor = root.querySelector<HTMLElement>('.ui-input-editor');
    if (!inputEditor) return;

    const highlightColor = '#e0f7fa';
    const defaultColor = '#ffffff';

    const setHighlightedStyle = () => inputEditor.style.backgroundColor = highlightColor;
    const resetStyle = () => inputEditor.style.backgroundColor = defaultColor;

    inputEditor.addEventListener('dragenter', setHighlightedStyle);
    inputEditor.addEventListener('dragover', (event) => {
      event.preventDefault();
      setHighlightedStyle();
    });

    inputEditor.addEventListener('dragleave', resetStyle);
    inputEditor.addEventListener('drop', (event) => {
      event.preventDefault();
      resetStyle();
    });
  }

  private static removeOptionsIcon(root: HTMLElement): void {
    const optionsIcon = root.querySelector<HTMLElement>(
      '.ui-input-options .grok-icon.fal.fa-cloud-upload'
    );
    optionsIcon?.remove();
  }
}

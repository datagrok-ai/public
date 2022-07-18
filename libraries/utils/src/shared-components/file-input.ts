import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subject} from 'rxjs';
import '../../css/shared-components.css';

export const EXCEL_BLOB_TYPE = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet';

export class FileInput {
  // events to emit
  public uploadedFile$ = new Subject<File | null>();

  public root = ui.div();

  private visibleInput = ui.stringInput('Input file', this.initialText, null);
  private hiddenInput = document.createElement('input');
  private icon = ui.iconFA('cloud-upload', () => this.hiddenInput.click());

  constructor(
    public initialText = 'Drag-n-drop here',
    public initialValue: File | null = null,
    public onValueChanged: Function | null = null,
    public fileType = EXCEL_BLOB_TYPE
  ) {
    this.draw();

    if (onValueChanged)
      this.uploadedFile$.subscribe((newValue: any) => onValueChanged(newValue));
  }

  private draw() {
    this.clear();

    const createFileArea = () => {
      (this.visibleInput.input as HTMLInputElement).readOnly = true;
      this.visibleInput.input.classList.add('default');

      const handleFiles = (files: FileList) => {
        if (files.length > 1) {
          this.visibleInput.input.classList.remove('success');
          this.visibleInput.input.classList.add('error');
          const newIcon = ui.iconFA('redo', () => this.hiddenInput.click());
          this.icon.replaceWith(newIcon);
          this.icon = newIcon;
          throw new Error('Please specify single input file');
        }

        this.visibleInput.value = files[0].name;
        if (files[0].type !== this.fileType) {
          this.visibleInput.input.classList.remove('success');
          this.visibleInput.input.classList.add('error');
          const newIcon = ui.iconFA('redo', () => this.hiddenInput.click());
          this.icon.replaceWith(newIcon);
          this.icon = newIcon;
          throw new Error('File type is not supported');
        }

        this.visibleInput.input.classList.add('success');
        const newIcon = ui.iconFA('times', () => this.reset());
        this.icon.replaceWith(newIcon);
        this.icon = newIcon;
        this.uploadedFile$.next(files[0]);
      };

      this.visibleInput.root.classList.add('drop-area');
      this.visibleInput.root.style.width = '100%';

      // hidden input to handle file dialog

      this.hiddenInput.type = 'file';
      this.hiddenInput.onchange = (e) => {
        //@ts-ignore
        const files: FileList = e.target.files;
        handleFiles(files);
      };

      // Prevent default drag behaviors
      ['dragenter', 'dragover', 'dragleave', 'drop'].forEach((eventName) => {
        this.visibleInput.input.addEventListener(eventName, (e) => {
          e.preventDefault();
          e.stopPropagation();
        }, false);
      });

      ['dragenter', 'dragover'].forEach((eventName) => {
        this.visibleInput.input.addEventListener(eventName, () => {
          this.visibleInput.input.classList.add('drag-n-dropping');
        }, false);
      });

      ['dragleave', 'drop'].forEach((eventName) => {
        this.visibleInput.input.addEventListener(eventName, () => {
          this.visibleInput.input.classList.remove('drag-n-dropping');
        }, false);
      });

      // Handle dropped files
      this.visibleInput.input.addEventListener('drop', async (e: DragEvent) => {
        const dt = e.dataTransfer;
        const files = dt!.files;
        handleFiles(files);
      }, false);

      // Pass clicks to hidden file input
      this.visibleInput.input.addEventListener('click', () => this.hiddenInput.click(), false);
      this.visibleInput.root.append(ui.div(this.icon, 'icon'));

      return this.visibleInput;
    };

    ['drag', 'dragenter'].forEach((eventName) => {
      document.getElementsByClassName('layout-workarea')[0].addEventListener(eventName, (e) => {
        e.stopPropagation();
      }, false);
    });

    this.root.style.display = 'flex';
    this.root.append(createFileArea().root);
  }

  private reset() {
    this.visibleInput.input.classList.value = 'ui-input-editor default';
    const newIcon = ui.iconFA('cloud-upload', () => this.hiddenInput.click());
    this.icon.replaceWith(newIcon);
    this.icon = newIcon;
    this.visibleInput.value = 'Drag-n-drop here';
    this.uploadedFile$.next(null);
  }

  private clear() {
    ui.empty(this.root);
  }
}

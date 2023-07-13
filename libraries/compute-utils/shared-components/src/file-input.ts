/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BehaviorSubject, Subject} from 'rxjs';
import Validation from './validation';
import '../css/file-input.css';
import {FuncCallInput, InputWrapper} from './FuncCallInput';

export const EXCEL_BLOB_TYPE = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet';

export class FileInput implements FuncCallInput<File>, InputWrapper<any> {
  // events to emit
  public uploadedFile$ = new Subject<File | null>();
  // legacy compatibility
  public onFileUploaded = this.uploadedFile$;

  public notify = true;
  public value: File | null = null;

  public get enabled() {
    return this.primaryInput.enabled;
  }

  public set enabled(val: boolean) {
    this.primaryInput.enabled = val;
  }

  public get input() {
    return this.primaryInput;
  }

  public get captionLabel() {
    return this.primaryInput.captionLabel;
  }

  // HTML root of component
  public root = ui.div();
  // Validation object
  public validation = new Validation<File>();

  public primaryInput = ui.stringInput('Input file', this.initialText, null);
  private hiddenInput = document.createElement('input');
  private icon = ui.iconFA('cloud-upload', () => this.hiddenInput.click());

  private isValid = new BehaviorSubject<boolean>(true);

  constructor(
    public initialText = 'Drag-n-drop here',
    public initialValue: File | null = null,
    public onValueChanged: Function | null = null,
    public fileType: string | null = EXCEL_BLOB_TYPE,
  ) {
    this.draw();

    this.isValid.subscribe((isValid) => {
      if (isValid) {
        this.primaryInput.input.classList.remove('error');
        this.primaryInput.input.classList.add('success');
        const newIcon = ui.iconFA('times', () => this.reset(), 'Reset uploaded file');
        this.icon.replaceWith(newIcon);
        this.icon = newIcon;
      } else {
        this.primaryInput.input.classList.remove('success');
        this.primaryInput.input.classList.add('error');
        const newIcon = ui.iconFA('redo', () => this.hiddenInput.click(), 'Re-upload a file');
        this.icon.replaceWith(newIcon);
        this.icon = newIcon;
      }
    });

    this.uploadedFile$.subscribe((newValue: File | null) => {
      if (this.onValueChanged)
        this.onValueChanged(newValue);

      this.value = newValue;
    });
  }

  public onInput(cb: Function) {
    return this.uploadedFile$.subscribe((file) => {
      if (this.notify)
        cb(file);
    });
  }

  private draw() {
    this.clear();

    const createFileArea = () => {
      (this.primaryInput.input as HTMLInputElement).readOnly = true;
      this.primaryInput.input.classList.add('default');

      const handleFiles = async (files: FileList | null) => {
        if (!files || !files.length || !this.enabled)
          return;

        if (files.length > 1) {
          this.isValid.next(false);
          throw new Error('Please specify single input file');
        }

        this.primaryInput.value = files[0].name;
        if (this.fileType && files[0].type !== this.fileType) {
          this.isValid.next(false);
          throw new Error('File type is not supported');
        }

        this.isValid.next(await this.validation.validate(files[0]));

        if (this.isValid.value)
          this.uploadedFile$.next(files[0]);
        else
          this.uploadedFile$.next(null);
      };

      this.primaryInput.root.classList.add('fi-drop-area');
      this.primaryInput.root.style.width = '100%';

      // hidden input to handle file dialog

      this.hiddenInput.type = 'file';
      this.hiddenInput.addEventListener('change', async () => await handleFiles(this.hiddenInput.files), false);

      // Prevent default drag behaviors
      ['dragenter', 'dragover', 'dragleave', 'drop'].forEach((eventName) => {
        this.primaryInput.input.addEventListener(eventName, (e) => {
          e.preventDefault();
          e.stopPropagation();
        }, false);
      });

      ['dragenter', 'dragover'].forEach((eventName) => {
        this.primaryInput.input.addEventListener(eventName, () => {
          this.primaryInput.input.classList.add('drag-n-dropping');
        }, false);
      });

      ['dragleave', 'drop'].forEach((eventName) => {
        this.primaryInput.input.addEventListener(eventName, () => {
          this.primaryInput.input.classList.remove('drag-n-dropping');
        }, false);
      });

      // Handle dropped files
      this.primaryInput.input.addEventListener('drop', async (e: DragEvent) => {
        const dt = e.dataTransfer;
        const files = dt!.files;
        handleFiles(files);
      }, false);

      // Pass clicks to hidden file input
      this.primaryInput.input.addEventListener('click', () => this.hiddenInput.click(), false);
      this.primaryInput.root.append(ui.div(this.icon, 'icon'));

      return this.primaryInput;
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
    this.primaryInput.input.classList.value = 'ui-input-editor default';
    const newIcon = ui.iconFA('cloud-upload', () => this.hiddenInput.click(), 'Choose a file to upload');
    this.icon.replaceWith(newIcon);
    this.icon = newIcon;
    this.primaryInput.value = 'Drag-n-drop here';
    this.hiddenInput.value = '';
    this.uploadedFile$.next(null);
  }

  private clear() {
    ui.empty(this.root);
  }
}

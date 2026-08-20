/* Private-key input: a key file goes in, its content is held masked (port of `rsa_input.dart`).
   dg, because a key dragged out of the platform's file browser is read with `grok.dapi.files`. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Input, InputOptions} from '../core/input-base.js';
import {button, span} from '../core/elements.js';
import {TextInput} from '../components/inputs/text-input.js';
import {buttonWithIcon} from '../components/actions/buttons.js';
import {dropZone} from './file-input.js';

export interface RsaInputOptions extends InputOptions<string | null> {
  /** Extensions the picker offers and a dropped file must carry; `.pem,.der` by default. */
  accept?: string;
}

/** Upload button and drop zone until a key is in, a masked field with a Change button afterwards.
 * The value is the key itself: the file's text when it opens with a PEM header, BASE64 of its
 * bytes when it does not (`rsa_input.dart:119-129`) — so the same input carries armoured and DER
 * keys. Files dropped from the platform's own file browser are read through
 * `grok.dapi.files.readAsText` and must be PEM, as in Dart (`:33-45`). */
export class RsaInput extends Input<string | null, RsaInputOptions> {
  /** What the platform shows in place of a stored secret (`credentials.dart:13`); js-api does not
   * surface the constant. */
  static readonly PLACEHOLDER = '_____________';
  private static readonly PEM = '-----BEGIN';

  private _picker!: HTMLInputElement;

  constructor(options: RsaInputOptions = {}) {
    super(options, null);
    this.root.dataset.u2 = 'rsa-input';
  }

  /** A key file's content as this input stores it: PEM text as it stands, anything else BASE64
   * (`rsa_input.dart:125-127`). */
  static content(bytes: Uint8Array): string {
    const text = new TextDecoder().decode(bytes).trim();
    return text.startsWith(RsaInput.PEM) ? text : RsaInput.base64(bytes);
  }

  static base64(bytes: Uint8Array): string {
    let binary = '';
    // chunked: String.fromCharCode(...bytes) overflows the argument stack on a real key file
    for (let at = 0; at < bytes.length; at += 0x8000)
      binary += String.fromCharCode.apply(null, [...bytes.subarray(at, at + 0x8000)]);
    return btoa(binary);
  }

  protected createEditor(): HTMLElement {
    const picker = document.createElement('input');
    this._picker = picker;
    picker.type = 'file';
    picker.accept = this._accept;
    picker.hidden = true;
    this._listen(picker, 'change', () => {
      const file = picker.files?.[0];
      if (file)
        this.accept(file);
    });

    const upload = document.createElement('div');
    upload.className = 'u2-rsa-upload';
    upload.append(picker, buttonWithIcon('Upload', 'paperclip', () => picker.click()),
      span(`or drop a ${this._accept} file`, 'u2-rsa-hint'));

    const masked = new TextInput({password: true, value: RsaInput.PLACEHOLDER, inline: true});
    masked.root.querySelector('input')!.readOnly = true;
    const maskedBox = document.createElement('div');
    maskedBox.className = 'u2-rsa-masked';
    maskedBox.append(masked.root, button('Change', () => this._change()));

    const editor = document.createElement('div');
    editor.className = 'u2-rsa';
    editor.append(upload, maskedBox);
    this.own(dropZone(editor, (files) => this._drop(files)));
    this._makeDroppable(editor);

    this.effect(() => {
      const stored = (this.value.value ?? '') !== '';
      upload.hidden = stored;
      maskedBox.hidden = !stored;
    });
    this.addValidator((v) => !this.nullable && (v ?? '') === '' ? 'Can\'t be empty.' : null);
    return editor;
  }

  /** Reads a local file and stores it; a file whose extension is not accepted is refused with a
   * message rather than silently ignored. */
  async accept(file: File): Promise<void> {
    if (!this._accepts(file.name)) {
      grok.shell.error(`Expected a ${this._accept} file.`);
      return;
    }
    const content = RsaInput.content(new Uint8Array(await file.arrayBuffer()));
    if (!this.scope.isDisposed)
      this.value.value = content;
  }

  private _drop(files: FileList): void {
    if (files.length > 1)
      grok.shell.error('Multiple files are not supported.');
    else
      this.accept(files[0]);
  }

  private get _accept(): string {
    return this.options.accept ?? '.pem,.der';
  }

  private _accepts(name: string): boolean {
    return this._accept.split(',').some((ext) => name.toLowerCase().endsWith(ext.trim().toLowerCase()));
  }

  private _change(): void {
    this.value.value = null;
    this._picker.value = '';
    this._picker.click();
  }

  /** The platform's drag channel, for keys dragged out of the file browser: those arrive as
   * `FileInfo`s, not as browser files, and only PEM text is taken (`rsa_input.dart:33-45`).
   * The registration is dart-side and has no counterpart to undo it, so it outlives dispose — a
   * drop that lands after it is read and thrown away rather than written. */
  private _makeDroppable(el: HTMLElement): void {
    ui.makeDroppable(el, {
      acceptDrop: (x) => this._isKeyFile(x),
      doDrop: async (payload) => {
        const dropped = RsaInput.dragged(payload);
        if (!this._isKeyFile(dropped))
          return;
        const content = await grok.dapi.files.readAsText(dropped).catch(() => null);
        if (this.scope.isDisposed)
          return;
        if (content === null)
          grok.shell.error('Could not read the key file.');
        else if (content.startsWith(RsaInput.PEM))
          this.value.value = content;
        else
          grok.shell.error('Invalid key file');
      },
    });
  }

  /** What the drop carries: js-api hands the callback a `DragDropArgs` (`events.ts:309-320`),
   * platform builds before it handed the dragged object itself. */
  static dragged(payload: any): any {
    return payload !== null && typeof payload === 'object' && 'dragObject' in payload ?
      payload.dragObject : payload;
  }

  /** A key file out of the platform's file browser — anything else on the drag channel is somebody
   * else's and is left alone. The connection is what makes it readable through `dapi`, and Dart
   * asks for it too (`rsa_input.dart:32`). */
  private _isKeyFile(x: any): x is DG.FileInfo {
    return x instanceof DG.FileInfo && x.isFile && RsaInput._connected(x) && this._accepts(x.fileName);
  }

  /** `FileInfo.connection` throws on every platform whose js-api bundle still passes the wrapped
   * handle to the accessor (`table-info.ts:73`), which would reject every drop. A connection that
   * cannot be read is taken as present and the read above decides. */
  private static _connected(x: DG.FileInfo): boolean {
    try {
      return x.connection != null;
    } catch {
      return true;
    }
  }

  private _listen(el: HTMLElement, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}

export function rsaInput(label: string, options: RsaInputOptions = {}): RsaInput {
  return new RsaInput({...options, label});
}

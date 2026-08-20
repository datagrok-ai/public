/* File pickers over the platform's file shares. dg, not core: a typed path is checked with
   `grok.dapi.files`, a remote one is namespaced by its connection, and the value is a platform
   `FileInfo` (port of `file_input.dart`). */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Input, InputOptions} from '../core/input-base.js';
import {signal, Signal} from '../core/signals.js';
import {AsyncSource, AsyncState} from '../core/async-source.js';
import {iconButton} from '../components/actions/buttons.js';

export type FileInputMode = 'file' | 'folder' | 'any';

export interface FileInputOptions extends InputOptions<DG.FileInfo | null> {
  /** What the path must resolve to; `file` by default (Dart `FileInputMode`). */
  mode?: FileInputMode;
  /** Share the path is relative to: its `nqName` prefixes the path before the check, and the
   * local-file affordance goes away (`file_input.dart:54,136-137`). */
  connection?: DG.DataConnection;
  placeholder?: string;
}

/** Path box with a local-file picker and a drop zone, validated against the file share as it is
 * typed. The platform's three validator states are kept (`file_input.dart:34-40`): empty, still
 * checking, and does-not-exist; the check itself runs through {@link AsyncSource}, so a keystroke
 * cancels the pending one and a late answer to an old path is dropped instead of overwriting a
 * newer verdict.
 *
 * Divergences from Dart, all deliberate: the message names what the mode asks for (`Folder does
 * not exist.` rather than `File …`) and never repeats the path, which the box shows anyway; the
 * root of a connection resolves as existing
 * with no `FileInfo` and leaves the value alone, since js-api has no way to build a directory
 * `FileInfo` (`:147-149` does it dart-side), where the platform input throws out of `firstWhere`
 * (`:156`); and an empty box is invalid only when the input is not nullable, where
 * `file_input.dart:22` defaults `nullable` to false and so always demands a path. */
export class FileInput extends Input<DG.FileInfo | null, FileInputOptions> {
  // all assigned in createEditor, which the base constructor runs before field initializers
  private _field!: HTMLInputElement;
  private _text!: Signal<string>;
  private _busy!: Signal<boolean>;
  private _problem!: Signal<string | null>;
  private _source!: AsyncSource<DG.FileInfo | null>;
  // bumped by everything that changes what the box holds; a check or a read whose turn has passed
  // is ignored rather than allowed to overwrite a newer answer
  private _turn!: number;
  private _queried!: number;
  // up while the input writes its own value, so `_show` knows the write is not a new subject
  private _echo!: boolean;

  constructor(options: FileInputOptions = {}) {
    super(options, null);
    this.root.dataset.u2 = 'file-input';
  }

  /** The path the platform is asked about: trailing slash dropped, and namespaced by the
   * connection unless it already is (`file_input.dart:134-137`). Null for an empty box. */
  static qualify(text: string, nqName?: string): string | null {
    if (text === '')
      return null;
    let path = text.endsWith('/') ? text.substring(0, text.lastIndexOf('/')) : text;
    if (nqName !== undefined && !path.startsWith(nqName))
      path = `${nqName}${path.startsWith('/') ? path : `/${path}`}`;
    return path;
  }

  protected createEditor(): HTMLElement {
    this._turn = 0;
    this._queried = 0;
    this._echo = false;
    this._text = signal('');
    this._busy = signal(false);
    this._problem = signal<string | null>(null);
    this._source = new AsyncSource<DG.FileInfo | null>((path) => this._resolve(path), {debounceMs: 250});
    this.own(() => this._source.dispose());
    // a local read still running when the input goes away has nothing left to write into
    this.own(() => this._turn++);

    const field = document.createElement('input');
    this._field = field;
    field.type = 'text';
    field.autocomplete = 'off';
    field.className = 'u2-file-path';
    field.placeholder = this.options.placeholder ?? '';
    this.effect(() => this._show(this.value.value));
    this._listen(field, 'input', () => {
      this._queried = ++this._turn;
      this._text.value = field.value;
      this._source.query(field.value);
    });

    const editor = document.createElement('div');
    editor.className = 'u2-file';
    editor.append(field);
    if (this.options.connection === undefined && this.options.mode !== 'folder')
      this._addPicker(editor);
    this.own(dropZone(editor, (files) => this._drop(files)));

    this.effect(() => this._onState(this._source.state.value));
    this.addValidator(() => this._check());
    return editor;
  }

  private get _noun(): string {
    return this.options.mode === 'folder' ? 'Folder' : 'File';
  }

  private _check(): string | null {
    const text = this._text.value;
    if (text === '')
      return this.nullable ? null : 'Can\'t be empty.';
    if (this._busy.value)
      return `${this._noun} ${text} is loading.`;
    return this._problem.value;
  }

  private _onState(state: AsyncState<DG.FileInfo | null>): void {
    if (state.kind === 'idle')
      return;
    if (state.kind === 'loading') {
      this._busy.value = true;
      return;
    }
    // a file picked or dropped since the check started is what the box holds now
    if (this._queried !== this._turn)
      return;
    this._busy.value = false;
    if (state.kind === 'error') {
      this._problem.value = state.message;
      return;
    }
    if (state.kind === 'empty') {
      this._problem.value = `${this._noun} does not exist.`;
      return;
    }
    this._problem.value = null;
    const found = state.items[0];
    if (found !== null)
      this._publish(found);
  }

  /** The path as the platform spells it — relative to the connection when there is one
   * (`file_input.dart:91-95`). A value written from outside is a new subject: the verdict on the
   * old path goes with it, and a check still in flight cannot re-assert one, as the platform
   * setter clears and revalidates (`file_input.dart:81-88`). */
  private _show(file: DG.FileInfo | null): void {
    const text = file === null ? '' :
      this.options.connection !== undefined ? (file.path ? file.path : '/') :
        file.fullPath ?? file.fileName;
    this._field.value = text;
    this._text.value = text;
    if (this._echo)
      return;
    this._turn++;
    this._busy.value = false;
    this._problem.value = null;
  }

  private _publish(file: DG.FileInfo | null): void {
    this._echo = true;
    try {
      this.value.value = file;
    } finally {
      this._echo = false;
    }
  }

  private async _resolve(text: string): Promise<(DG.FileInfo | null)[]> {
    const path = FileInput.qualify(text, this.options.connection?.nqName);
    if (path === null || !await grok.dapi.files.exists(path))
      return [];
    const found = await this._find(path);
    if (found === null)
      return this._isConnectionRoot(path) ? [null] : [];
    const mode = this.options.mode ?? 'file';
    if (mode === 'file' && !found.isFile || mode === 'folder' && !found.isDirectory)
      return [];
    return [found];
  }

  /** The one path that resolves with no `FileInfo` to show for it (`file_input.dart:141-143`).
   * Anything else the listing does not carry did not resolve, and says so — a share that answers
   * every call with an empty list would otherwise leave the input silently valid. */
  private _isConnectionRoot(path: string): boolean {
    const nqName = this.options.connection?.nqName;
    return nqName !== undefined && path.length === nqName.length + 1;
  }

  /** js-api has no per-path lookup, so the entry comes out of its directory listing — which is
   * where the platform input gets the `FileInfo` it stores too (`file_input.dart:151-156`). The
   * extension filter is for file mode alone: a directory whose own name carries a dot (`geo/v1.2`)
   * would otherwise be looked for among `.2` files. */
  private async _find(path: string): Promise<DG.FileInfo | null> {
    const slash = path.lastIndexOf('/');
    const dot = path.lastIndexOf('.');
    const ext = (this.options.mode ?? 'file') === 'file' && dot > slash ?
      path.substring(dot + 1) : null;
    const list = await grok.dapi.files.list(slash < 0 ? '' : path.substring(0, slash + 1), false, ext);
    return list.find((f) => f.fullPath === path) ?? null;
  }

  /** The local-file affordance on the shared options rail, where the platform input puts it too
   * (`file_input.dart:54-55` through `InputBase.addOptions`); the file element that does the
   * picking is hidden and stays out of the rail, which only carries what shows. */
  private _addPicker(editor: HTMLElement): void {
    const picker = document.createElement('input');
    picker.type = 'file';
    picker.hidden = true;
    this._listen(picker, 'change', () => {
      const file = picker.files?.[0];
      if (file)
        this._accept(file);
    });
    editor.append(picker);
    this.addOptions(iconButton('folder-open', () => picker.click(), {tooltip: 'Open file'}));
  }

  private _drop(files: FileList): void {
    if (files.length > 1)
      grok.shell.error('Multiple files are not supported.');
    else
      this._accept(files[0]);
  }

  /** A local file becomes an in-memory `FileInfo`, as the platform's reader does
   * (`file_input.dart:108-128`). A read still running when the box is typed into is abandoned. */
  private async _accept(file: File): Promise<void> {
    const turn = ++this._turn;
    this._busy.value = true;
    this._field.value = file.name;
    this._text.value = file.name;
    try {
      const bytes = new Uint8Array(await file.arrayBuffer());
      if (turn !== this._turn)
        return;
      this._problem.value = null;
      this._publish(DG.FileInfo.fromBytes(file.name, bytes));
    } catch (e) {
      if (turn === this._turn)
        this._problem.value = `${this._noun} ${file.name} could not be read.`;
    } finally {
      if (turn === this._turn)
        this._busy.value = false;
    }
  }

  private _listen(el: HTMLElement, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}

export function fileInput(label: string, options: FileInputOptions = {}): FileInput {
  return new FileInput({...options, label});
}

/** Files dropped on `el`, with the browser's own navigate-to-the-file default and the drag
 * chrome handled; returns the cleanup its caller owns. */
export function dropZone(el: HTMLElement, onFiles: (files: FileList) => void): () => void {
  const mark = (e: Event, dropping: boolean) => {
    e.preventDefault();
    e.stopPropagation();
    el.classList.toggle('u2-file-dropping', dropping);
  };
  const handlers: {[type: string]: (e: Event) => void} = {
    dragenter: (e) => mark(e, true),
    dragover: (e) => mark(e, true),
    dragleave: (e) => mark(e, false),
    drop: (e) => {
      mark(e, false);
      const files = (e as DragEvent).dataTransfer?.files;
      if (files && files.length > 0)
        onFiles(files);
    },
  };
  for (const type in handlers)
    el.addEventListener(type, handlers[type]);
  return () => {
    for (const type in handlers)
      el.removeEventListener(type, handlers[type]);
  };
}

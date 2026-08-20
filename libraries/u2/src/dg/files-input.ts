/* Multi-file input: rows with per-file progress while the bytes arrive, one summary line when
   they have (port of `files_input.dart`). dg, because every loaded file becomes a platform
   `FileInfo`. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Input, InputOptions} from '../core/input-base.js';
import {signal, Signal} from '../core/signals.js';
import {ProgressBar} from '../components/progress.js';
import {iconButton} from '../components/buttons.js';
import {dropZone} from './file-input.js';

export interface FilesInputOptions extends InputOptions<DG.FileInfo[]> {
  /** Extensions the picker offers and dropped files must carry, e.g. `['.csv', '.sdf']`. */
  accept?: string[];
}

interface FileLoad {
  key: string;
  row: HTMLElement;
  status: HTMLElement;
  info: DG.FileInfo | null;
  error: string | null;
  loading: boolean;
  cancelled: boolean;
  cancel: () => void;
}

/** Picks or accepts a drop of several files at once: each gets a row with its own percentage and
 * a ✕ that cancels the read or drops the file, and the value — the loaded `FileInfo`s — appears
 * once every read has settled (`files_input.dart:250`). Files are keyed by name, so re-adding one
 * replaces it rather than duplicating it (`:100-102`).
 *
 * Two deliberate divergences from Dart: the value is written whenever the reads settle, failures
 * and all (the platform input withholds it while invalid, which hides a partial load from a form
 * that only reads the value), and an empty selection is invalid only when the input is not
 * nullable — where `files_input.dart:37-38` requires a file regardless. The bytes are read off
 * `File.stream()` rather than a `FileReader`, which is what makes cancel and progress work the
 * same in a browser and in a test. */
export class FilesInput extends Input<DG.FileInfo[], FilesInputOptions> {
  // all assigned in createEditor, which the base constructor runs before field initializers
  private _items!: Signal<FileLoad[]>;
  private _list!: HTMLElement;
  private _summary!: HTMLElement;
  private _progress!: ProgressBar;
  private _echo!: boolean;

  constructor(options: FilesInputOptions = {}) {
    super(options, []);
    this.root.dataset.u2 = 'files-input';
  }

  /** Adds files as if they had been dropped on the input — the picker and the drop zone path. */
  add(files: ArrayLike<File>): void {
    const items = [...this._items.peek()];
    for (let i = 0; i < files.length; i++) {
      const file = files[i];
      if (!this._accepts(file.name)) {
        grok.shell.warning(`${file.name} is not ${this.options.accept!.join(' or ')}.`);
        continue;
      }
      const existing = items.findIndex((x) => x.key === file.name);
      if (existing >= 0)
        this._drop(items[existing], items, existing);
      const item = this._item(file.name);
      items.push(item);
      this._load(file, item);
    }
    this._items.value = items;
  }

  protected createEditor(): HTMLElement {
    this._echo = false;
    this._items = signal<FileLoad[]>([]);
    this._progress = new ProgressBar({showPercent: true});
    this._progress.root.hidden = true;

    this._summary = document.createElement('span');
    this._summary.className = 'u2-files-summary';
    this._summary.hidden = true;
    this._list = document.createElement('div');
    this._list.className = 'u2-files-list';

    const editor = document.createElement('div');
    editor.className = 'u2-files';
    editor.append(this._summary, this._progress.root, this._list);
    this._addPicker();
    this.own(dropZone(editor, (files) => this.add(files)));
    // a read still running when the input goes away would settle onto a disposed component
    this.own(() => {
      for (const item of this._items.peek())
        this._cancel(item);
    });

    this.effect(() => this._render(this._items.value));
    this.effect(() => this._adopt(this.value.value));
    this.addValidator(() => this._check());
    return editor;
  }

  private _check(): string | null {
    const items = this._items.value;
    if (items.some((i) => i.loading))
      return 'Files are still loading.';
    if (items.some((i) => i.error !== null))
      return 'Some files failed to load.';
    if (items.length === 0 && !this.nullable)
      return 'At least one file is required.';
    return null;
  }

  private _accepts(name: string): boolean {
    const accept = this.options.accept;
    return accept === undefined || accept.some((ext) => name.toLowerCase().endsWith(ext.toLowerCase()));
  }

  private _render(items: FileLoad[]): void {
    this._list.replaceChildren(...items.map((i) => i.row));
    const loading = items.filter((i) => i.loading).length;
    const failed = items.filter((i) => i.error !== null).length;
    const loaded = items.length - loading - failed;
    this._progress.root.hidden = loading === 0;
    this._summary.hidden = loading > 0 || items.length === 0;
    this._summary.textContent = failed > 0 ? `${loaded} loaded, ${failed} failed` : `${loaded} loaded`;
    if (loading > 0)
      this._progress.value.value = loaded / items.length;
    this.refreshEnabled();
  }

  /** A value written from outside replaces the rows, as the platform input rebuilds its tags
   * (`files_input.dart:194-207`); the rows written by a finished read do not come back through
   * here. */
  private _adopt(files: DG.FileInfo[]): void {
    if (this._echo)
      return;
    for (const item of this._items.peek())
      this._cancel(item);
    this._items.value = files.map((f) => {
      const item = this._item(f.fileName ?? f.path);
      item.info = f;
      item.loading = false;
      item.status.textContent = '';
      return item;
    });
  }

  private _item(name: string): FileLoad {
    const row = document.createElement('div');
    row.className = 'u2-files-row';
    const label = document.createElement('span');
    label.className = 'u2-files-name';
    label.textContent = name;
    label.title = name;
    const status = document.createElement('span');
    status.className = 'u2-files-status';
    status.textContent = '0%';
    const item: FileLoad = {key: name, row, status, info: null, error: null, loading: true,
      cancelled: false, cancel: () => {}};
    row.append(label, status, iconButton('times', () => this._remove(item), {tooltip: 'Remove'}));
    return item;
  }

  private async _load(file: File, item: FileLoad): Promise<void> {
    const reader = file.stream().getReader();
    // cancelling a reader whose stream already failed rejects, and nothing awaits it here
    item.cancel = () => void reader.cancel().catch(() => {});
    const chunks: Uint8Array[] = [];
    let loaded = 0;
    try {
      for (;;) {
        const {done, value} = await reader.read();
        if (item.cancelled)
          return;
        if (done || value === undefined)
          break;
        chunks.push(value);
        loaded += value.length;
        item.status.textContent = file.size > 0 ? `${Math.round(loaded / file.size * 100)}%` : '';
      }
      item.info = DG.FileInfo.fromBytes(file.name, FilesInput._join(chunks, loaded));
    } catch (e) {
      item.error = e instanceof Error ? e.message : String(e);
      grok.shell.error(`Failed to load ${file.name}`);
    }
    item.loading = false;
    item.status.textContent = item.error ?? '';
    item.row.classList.toggle('u2-files-failed', item.error !== null);
    item.row.title = item.error ?? '';
    this._settle();
  }

  private static _join(chunks: Uint8Array[], length: number): Uint8Array {
    const bytes = new Uint8Array(length);
    let at = 0;
    for (const chunk of chunks) {
      bytes.set(chunk, at);
      at += chunk.length;
    }
    return bytes;
  }

  private _remove(item: FileLoad): void {
    const items = [...this._items.peek()];
    const at = items.indexOf(item);
    if (at < 0)
      return;
    this._drop(item, items, at);
    this._items.value = items;
    this._settle();
  }

  private _drop(item: FileLoad, items: FileLoad[], at: number): void {
    this._cancel(item);
    items.splice(at, 1);
  }

  private _cancel(item: FileLoad): void {
    item.cancelled = true;
    item.cancel();
  }

  /** Re-runs the summary and the validators, and publishes the value once nothing is still
   * loading — the rows are the source of truth, the value follows them. */
  private _settle(): void {
    const items = this._items.peek();
    this._items.value = [...items];
    if (items.some((i) => i.loading))
      return;
    this._echo = true;
    try {
      this.value.value = items.map((i) => i.info).filter((f): f is DG.FileInfo => f !== null);
    } finally {
      this._echo = false;
    }
  }

  /** The picker on the shared options rail, where the platform input puts its own summary
   * affordance (`files_input.dart:44` through `InputBase.addOptions`). */
  private _addPicker(): void {
    const picker = document.createElement('input');
    picker.type = 'file';
    picker.multiple = true;
    picker.hidden = true;
    const accept = this.options.accept;
    if (accept !== undefined)
      picker.accept = accept.join(',');
    const onChange = () => {
      if (picker.files)
        this.add(picker.files);
      picker.value = '';
    };
    picker.addEventListener('change', onChange);
    this.own(() => picker.removeEventListener('change', onChange));
    this.addOptions(picker);
    this.addOptions(iconButton('folder-open', () => picker.click(), {tooltip: 'Open files'}));
  }
}

export function filesInput(label: string, options: FilesInputOptions = {}): FilesInput {
  return new FilesInput({...options, label});
}

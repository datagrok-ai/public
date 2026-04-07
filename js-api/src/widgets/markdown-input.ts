/**
 * MarkdownInput class (separate file for Quill lazy loading).
 * @module widgets/markdown-input
 */

import {Utils} from '../utils';
import {JsInputBase} from "./inputs-base";
import {MarkdownConfig} from "./types";

declare let DG: any;
declare let ui: any;


export class MarkdownInput extends JsInputBase<string> {
  // Quill.js instance - loaded from the core .min.js file
  private editor: any;
  private readonly _editorRoot: HTMLDivElement = ui.div([], 'markdown-input');
  private _textToSet?: string | null;

  private constructor(caption?: string) {
    super();
    this.root.classList.add('markdown-input-root');
    this.addCaption(caption ?? '');
  }

  static create(caption?: string, options?: MarkdownConfig): MarkdownInput {
    const input = new MarkdownInput(caption);
    ui.setUpdateIndicator(input.root, true, 'Loading markdown input...');
    Utils.loadJsCss([
      '/js/common/quill/quill.min.js',
      '/js/common/quill/quill.snow.css',
      '/js/common/quill/quilljs-markdown.min.js',
    ]).then(() => {
      //@ts-ignore
      input.editor = new Quill(input._editorRoot, {
        modules: {
          toolbar: [
            [{ header: [1, 2, false] }],
            ['bold', 'italic', 'underline', 'strike'],
            [{ list: 'ordered' }, { list: 'bullet' }],
            ['blockquote', 'code-block', 'image'],
          ],
        },
        theme: 'snow', // or 'bubble'
      });
      if (input._textToSet != null || (options && options.value != null && options.value !== '')) {
        //@ts-ignore
        input.editor.setText(input._textToSet ?? options.value);
        input._textToSet = null;
      }
      input.editor.on('text-change', (_: any, __: any, source: string) => {
        if (source === 'api')
          input.fireChanged();
        else if (source === 'user')
          input.fireInput();
      });
      //@ts-ignore
      new QuillMarkdown(input.editor, {
        syntax: true, // enables code blocks, etc.
        preview: true,
      });
      ui.setUpdateIndicator(input.root, false);
    }).catch((_) => {
      ui.setUpdateIndicator(input.root, false);
    });

    return input;
  }

  getInput(): HTMLElement { return this._editorRoot; }

  get inputType(): string { return 'markdown'; }

  get dataType(): string { return `${DG.TYPE.STRING}`; }

  get attachedImages(): string[] {
    const ops: { [key:string]: any; }[] = this.editor?.getContents()['ops'] ?? [];
    return ops
        .filter((op) => op['insert']?.constructor === Object
            && (op['insert']['image']?.startsWith('data:image') ?? false))
        .map((op) => op['insert']['image']);
  }

  getStringValue(): string {
    return this.editor?.getText() ?? '';
  }

  getValue(): string {
    return this.editor?.root?.innerHTML;
  }

  setStringValue(value: string): void {
    this.editor ? this.editor.setText(value): this._textToSet = value;
  }

  setValue(value: any): void {
    this.editor ? this.editor.setText(value): this._textToSet = value;
  }
}

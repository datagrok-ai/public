/**
 * CodeEditor class.
 * @module widgets/code-editor
 */

import {toJs} from "../wrappers";
import {Observable} from "rxjs";
import {observeStream} from "../events";
import {IDartApi} from "../api/grok_api.g";
import {ICodeEditorOptions} from "./types";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


/** Class for code input editor. */
export class CodeEditor {
  dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  static create(script = '', mode = 'javascript', placeholder = '', options?: ICodeEditorOptions): CodeEditor {
    return toJs(api.grok_CodeEditor(script, mode, placeholder, options?.root));
  }

  append(text: string): void {
    api.grok_CodeEditor_Append(this.dart, text);
  }

  async setReadOnly(value: boolean): Promise<void> {
    await api.grok_CodeEditor_SetReadOnly(this.dart, value);
  }

  get root(): HTMLElement { return api.grok_CodeEditor_Get_Root(this.dart); }

  get value(): string { return api.grok_CodeEditor_Get_Value(this.dart); }
  set value(x: string) { api.grok_CodeEditor_Set_Value(this.dart, x); }

  get onValueChanged(): Observable<any> { return observeStream(api.grok_CodeEditor_OnValueChanged(this.dart)); }
}

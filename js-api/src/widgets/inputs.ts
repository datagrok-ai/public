/**
 * Input classes: DateInput, ChoiceInput, TypeAhead, CodeInput.
 * @module widgets/inputs
 */

import {toDart, toJs} from "../wrappers";
import {Observable} from "rxjs";
import {IDartApi} from "../api/grok_api.g";
import dayjs from "dayjs";
import typeahead from 'typeahead-standalone';
import {Dictionary, typeaheadConfig} from 'typeahead-standalone/dist/types';
import {InputBase} from "./inputs-base";
import {CodeEditor} from "./code-editor";
import {TypeAheadConfig, CodeConfig} from "./types";

import '../../css/typeahead-input.css';

declare let ui: any;
const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


export class DateInput extends InputBase<dayjs.Dayjs | null> {
  declare dart: any;

  constructor(dart: any, onChanged: any = null) {
    super(dart, onChanged);
  }

  get value(): dayjs.Dayjs | null {
    const date = api.grok_DateInput_Get_Value(this.dart);
    return date == null ? date : dayjs(date);
  }
  set value(x: dayjs.Dayjs | null) { toDart(api.grok_DateInput_Set_Value(this.dart, x?.valueOf())); }
}

export class ChoiceInput<T> extends InputBase<T> {
  declare dart: any;

  constructor(dart: any, onChanged: any = null) {
    super(dart, onChanged);
  }

  get items(): T[] { return toJs(api.grok_ChoiceInput_Get_Items(this.dart)); }
  set items(s: T[]) { api.grok_ChoiceInput_Set_Items(this.dart, toDart(s)); }
}


export class TypeAhead extends InputBase {
  constructor(name: string, config: TypeAheadConfig) {
    const inputElement = ui.input.string(name, {value: ''});
    super(inputElement.dart);

    const typeAheadConfig: typeaheadConfig<Dictionary> = Object.assign(
      {input: <HTMLInputElement> this.input}, config);

    typeahead(typeAheadConfig);
    this._changeStyles();
  }

  _changeStyles() {
    this.root.getElementsByClassName('tt-list')[0].className = 'ui-input-list';
    this.root.getElementsByClassName('ui-input-list')[0].removeAttribute('style');
    this.root.getElementsByClassName('tt-input')[0].className = 'ui-input-editor';
    this.root.getElementsByClassName('typeahead-standalone')[0].classList.add('ui-input-root');
  }
}


export class CodeInput extends InputBase {
  declare dart: any;
  editor: CodeEditor;

  constructor(name: string, config?: CodeConfig) {
    const inputElement = ui.input.textArea(name);
    super(inputElement.dart);

    this.editor = CodeEditor.create(config?.script ?? '', config?.mode ?? 'javascript', config?.placeholder ?? '');
    this.input.style.display = 'none';
    this.root.classList.add('ui-input-code');
    this.root.append(this.editor.root);
  }

  get value(): string { return this.editor.value; }
  set value(x: string) { this.editor.value = x; }

  get onValueChanged(): Observable<any> { return this.editor.onValueChanged; }
}

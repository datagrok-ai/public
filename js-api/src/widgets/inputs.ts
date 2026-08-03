/**
 * Input classes: DateInput, ChoiceInput, MultiChoiceInput, TypeAhead, CodeInput.
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
import {TypeAheadConfig, TypeAheadCallbackSource, CodeConfig} from "./types";

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


export class MultiChoiceInput<T> extends InputBase<T[] | null> {
  declare dart: any;

  constructor(dart: any, onChanged: any = null) {
    super(dart, onChanged);
  }

  get items(): T[] { return toJs(api.grok_MultiChoiceInput_Get_Items(this.dart)); }
  set items(s: T[]) { api.grok_MultiChoiceInput_Set_Items(this.dart, toDart(s)); }
}


export class TypeAhead extends InputBase {
  constructor(name: string, config: TypeAheadConfig) {
    const inputElement = ui.input.string(name, {value: ''});
    super(inputElement.dart);

    if (typeof config.source === 'function') {
      this._initCallbackSource(config.source, config);
      return;
    }

    const typeAheadConfig: typeaheadConfig<Dictionary> = Object.assign(
      {input: <HTMLInputElement> this.input}, config) as typeaheadConfig<Dictionary>;

    typeahead(typeAheadConfig);
    this._changeStyles();
  }

  /** Self-contained dropdown for callback sources: typeahead-standalone only supports
   * local/prefetch/remote URL sources, so async callbacks (server-backed suggestions) get a
   * minimal debounced dropdown reusing the same CSS classes as the standalone path. */
  private _initCallbackSource(source: TypeAheadCallbackSource, config: TypeAheadConfig): void {
    const input = <HTMLInputElement> this.input;
    const wrapper = document.createElement('div');
    wrapper.className = 'typeahead-standalone ui-input-root';
    input.parentElement!.insertBefore(wrapper, input);
    wrapper.appendChild(input);
    const list = document.createElement('div');
    list.className = 'ui-input-list tt-hide';
    wrapper.appendChild(list);

    const minLength = config.minLength ?? 1;
    const limit = config.limit ?? 5;
    let timer: any = null;
    let requestId = 0;
    let items: Dictionary[] = [];
    let selectedIdx = -1;

    const hide = () => { list.classList.add('tt-hide'); selectedIdx = -1; };
    const pick = (ev: Event, i: number) => {
      const item = items[i];
      this.value = `${item.label ?? ''}`;
      hide();
      config.onSubmit?.(ev, item);
    };
    const render = () => {
      list.textContent = '';
      items.forEach((item, i) => {
        const el = document.createElement('div');
        el.className = i === selectedIdx ? 'tt-suggestion tt-selected' : 'tt-suggestion';
        el.textContent = `${item.label ?? ''}`;
        el.addEventListener('mousedown', (ev) => { ev.preventDefault(); pick(ev, i); });
        list.appendChild(el);
      });
      list.classList.toggle('tt-hide', items.length === 0);
    };

    input.addEventListener('input', () => {
      if (timer != null)
        clearTimeout(timer);
      const query = input.value;
      if (query.length < minLength) {
        items = [];
        hide();
        return;
      }
      timer = setTimeout(async () => {
        const id = ++requestId;
        let results: (string | Dictionary)[];
        try {
          results = await source(query);
        }
        catch (_) {
          results = [];
        }
        if (id !== requestId)
          return;
        items = results.slice(0, limit).map((x) => typeof x === 'string' ? {label: x} : x);
        selectedIdx = -1;
        render();
      }, config.debounceRemote ?? 100);
    });
    input.addEventListener('keydown', (ev) => {
      if (list.classList.contains('tt-hide'))
        return;
      if (ev.key === 'ArrowDown' || ev.key === 'ArrowUp') {
        selectedIdx = (selectedIdx + (ev.key === 'ArrowDown' ? 1 : items.length - 1)) % items.length;
        render();
        ev.preventDefault();
      }
      else if (ev.key === 'Enter' && selectedIdx >= 0) {
        pick(ev, selectedIdx);
        ev.preventDefault();
        ev.stopPropagation();
      }
      else if (ev.key === 'Escape')
        hide();
    });
    input.addEventListener('blur', () => hide());
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

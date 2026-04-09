/**
 * Dialog and InputForm classes.
 * @module widgets/forms
 */

import {toDart, toJs} from "../wrappers";
import {Observable} from "rxjs";
import {__obs, EventData, InputArgs, observeStream} from "../events";
import {Completer} from '../utils';
import {FuncCall} from "../functions";
import {IDartApi} from "../api/grok_api.g";
import {DartWidget, DartWrapper, Widget} from "./base";
import {InputBase} from "./inputs-base";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


/**
 * A non-modal dialog.
 * Sample: https://public.datagrok.ai/js/samples/ui/dialogs
 *
 * @example
 * ui.dialog('Windows')
 *   .add(ui.)
 *   .add(ui.span(['People of Earth, your attention, please… ']))
 *   .onOK(() => { grok.shell.info('OK!'); })
 *   .show();
 * */
  export class Dialog<Inputs extends Record<string, InputBase<any> > = {} > extends DartWidget {

  constructor(dart: any) {
    super(dart);
  }

  /** Creates a new dialog with the specified options. */
  static create(options?: { title?: string, helpUrl?: string, showHeader?: boolean, showFooter?: boolean } | string): Dialog {
    if (typeof options === 'string')
      return Dialog.create({ title: options });
    else
      return new Dialog(api.grok_Dialog(options?.title ?? '', options?.helpUrl, options?.showHeader ?? true, options?.showFooter ?? true));
  }

  /** When provided, adds a "?" icon to the dialog header on the right. */
  get helpUrl(): string { return api.grok_Dialog_Get_HelpUrl(this.dart); };
  set helpUrl(url: string) { api.grok_Dialog_Set_HelpUrl(this.dart, url); };

  /** Returns the title of a dialog. */
  get title(): string { return api.grok_Dialog_Get_Title(this.dart); };
  set title(t: string) { api.grok_Dialog_Set_Title(this.dart, t); };

  /** Returns a list of the dialog's inputs. */
  get inputs(): InputBase[] { return api.grok_Dialog_Get_Inputs(this.dart); }

  /** Returns an input with the specified caption, or throws an exception. */
  input(caption: string): InputBase {
    const input = this.inputs.find((f) => f.caption == caption);
    if (!input)
      throw `Input "${caption}" not found.`;
    return input;
  }

  /**
   * Sets the OK button handler, and shows the OK button
   * @param {Function} handler
   * @param {Object} options
   * @returns {Dialog} */
  onOK(handler: Function, options?: {closeOnEnter?: boolean}): Dialog<Inputs> {
    api.grok_Dialog_OnOK(this.dart, handler, options?.closeOnEnter ?? true);
    return this;
  }

  /**
   * Sets the OK button handler and returns a promise of the handler callback.
   * @param {Function} handler
   * @returns {Promise} */
  async awaitOnOK<T = any>(handler: () => Promise<T>): Promise<T> {
    let completer = new Completer<T>();
    this.onOK(() => {
      handler().then((res) => completer.complete(res))
               .catch((error) => completer.reject(error));
    });
    this.onCancel(() => {
      completer.reject();
    });
    return completer.promise;
  }

  /**
   * Sets the CANCEL button handler
   * @param {Function} handler
   * @returns {Dialog} */
  onCancel(handler: Function): Dialog<Inputs> {
    api.grok_Dialog_OnCancel(this.dart, handler);
    return this;
  }

  /** @returns {Observable} */
  get onClose(): Observable<any> {
    return __obs('d4-dialog-closed', this.dart);
  }

  // Using __obs is a recommended method. The below are obsolete and shall not be used:
  // onClose(handler) { api.grok_Dialog_OnClose(this.dart, handler); return this; }
  // onClose(handler) { let s = _sub(api.grok_Dialog_OnClose(this.dart, () => { handler(); s.cancel(); })); return this; }

  /** @returns {Dialog}
   * @param {{modal: boolean, fullScreen: boolean, center: boolean, centerAt: Element, x: number, y: number, width: number, height: number}|{}} options
   * */
  show(options?: { modal?: boolean; resizable?: boolean; fullScreen?: boolean; center?: boolean; centerAt?: Element; x?: number; y?: number; width?: number; height?: number; backgroundColor?: string; showNextTo?: HTMLElement}): Dialog<Inputs> {
    api.grok_Dialog_Show(this.dart, options?.modal, options?.resizable, options?.fullScreen, options?.center, options?.centerAt, options?.x, options?.y, options?.width, options?.height, options?.backgroundColor, options?.showNextTo);
    return this;
  }

  /** @returns {Dialog}
   * @param {boolean} fullScreen  */
  showModal(fullScreen: boolean): Dialog<Inputs> {
    api.grok_Dialog_Show(this.dart, true, null, fullScreen, false, null, null, null, null, null, null, null);
    return this;
  }

  /** Adds content to the dialog. using addInput() for inputs is preferred, as it provides better type safety.
   * @param {HTMLElement | Widget | InputBase} content
   * @returns {Dialog} */
  add(content: HTMLElement | Widget | InputBase): Dialog<Inputs> {
    api.grok_Dialog_Add(this.dart, toDart(content));
    return this;
  }

  /** Adds named input, which than can be used in combo with namedInputs, to get inputs in type-safe way*/
  addInput<K extends string, V extends InputBase<any>>(caption: K, input: V): Dialog<Record<K, V> & Inputs > {
    input.caption = caption;
    api.grok_Dialog_Add(this.dart, toDart(input));
    return this as unknown as Dialog<Record<K, V> & Inputs>;
  }

  /** gets typed inputs with captions as keys.*/
  get namedInputs() : Inputs {
    return this.inputs.reduce((acc, input) => {
      if (input.caption) {
        acc[input.caption] = input;
      }
      return acc;
    }, {} as Record<string, InputBase>) as Inputs;
  }

  /** Closes the dialog. */
  close(): void {
    api.grok_Dialog_Close(this.dart);
  }

  /** Returns command button with the specified text.
   * @param {string} text
   * @returns {HTMLButtonElement}
   * */
  getButton(text: string): HTMLButtonElement {
    return api.grok_Dialog_GetButton(this.dart, text);
  }

  /** Adds command button with the specified text.
   * @param {string} text
   * @param {Function} action
   * @param index
   * @param tooltip
   * @returns {Dialog}
   * */
  addButton(text: string, action: Function, index: number = 0, tooltip: any = null): Dialog<Inputs> {
    api.grok_Dialog_AddButton(this.dart, text, action, index, tooltip);
    return this;
  }

  /** Adds context action with the specified text.
   * @param {string} text
   * @param {Function} action
   * @returns {Dialog}
   * */
  addContextAction(text: string, action: Function): Dialog<Inputs> {
    api.grok_Dialog_AddContextAction(this.dart, text, action);
    return this;
  }

  /** Initializes the 'history' feature.
   * @param {Function} getInput - collects the input from UI into JSON-serializable object
   * @param {Function} applyInput - refreshes the UI according to input
   * */
  history(getInput: () => any, applyInput: (x: any) => void): void {
    api.grok_Dialog_History(this.dart, getInput, applyInput);
  }

  /** Initializes default history. */
  initDefaultHistory(): Dialog<Inputs> {
    api.grok_Dialog_InitDefaultHistory(this.dart);
    return this;
  }

  /** Clears the content. */
  clear() {
    api.grok_Dialog_Clear(this.dart);
  }

  /** Returns currently open dialogs. */
  static getOpenDialogs(): Dialog[] {
    return api.grok_Dialog_GetOpenDialogs();
  }

  /** Initializes the dialog properties from local storage. */
  initFromLocalStorage(): void {
    api.grok_Dialog_InitFromLocalStorage(this.dart);
  }
}


/** A form with multiple inputs inside */
export class InputForm extends DartWrapper {
  constructor(dart: any) { super(dart); }

  /** Creates an InputForm for the specified function call. */
  static async forFuncCall(funcCall: FuncCall, options?: { twoWayBinding?: boolean, skipDefaultInit?: boolean }): Promise<InputForm> {
    return new InputForm(await api.grok_InputForm_ForFuncCallAsync(funcCall.dart, options?.twoWayBinding ?? true, options?.skipDefaultInit ?? false));
  }

  static forInputs(inputs: InputBase[]): InputForm {
    inputs = inputs.filter((input) => input != null);
    return new InputForm(api.grok_InputForm_ForInputs(inputs.map((input) => input.dart)));
  }

  get root(): HTMLElement { return api.grok_InputForm_Get_Root(this.dart); };

  getInput(propertyName: string): InputBase { return toJs(api.grok_InputForm_GetInput(this.dart, propertyName)); }

  /** All inputs added to the form */
  get inputs(): InputBase[] {
    return api.grok_InputForm_GetInputs(this.dart);
  }

  get source(): any { return toJs(api.grok_InputForm_Get_Source(this.dart)); };

  set source(source: any) { api.grok_InputForm_Set_Source(this.dart, toDart(source)); };

  /** Occurs when user changes any input value in a form. */
  get onInputChanged(): Observable<EventData<InputArgs>> { return observeStream(api.grok_InputForm_OnInputChanged(this.dart)); }

  /** Occurs after the form is validated, no matter whether it is valid or not. */
  get onValidationCompleted(): Observable<any> { return observeStream(api.grok_InputForm_OnValidationCompleted(this.dart)); }

  /** Returns true if all inputs are valid. */
  get isValid(): boolean { return api.grok_InputForm_Get_IsValid(this.dart); }
}

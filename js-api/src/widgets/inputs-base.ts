/**
 * InputBase and JsInputBase classes.
 * @module widgets/inputs-base
 */

import {toDart, toJs} from "../wrappers";
import {Observable} from "rxjs";
import {Property} from "../entities";
import {Column} from "../dataframe";
import {IDartApi} from "../api/grok_api.g";
import * as d4 from '../api/d4.api.g';

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

declare let DG: any;


/** Input control base. Could be used for editing {@link Property} values as well.
 * The root is a div that consists of {@link captionLabel} and {@link input}.
 * */
export class InputBase<T = any> {
  dart: any;

  constructor(dart: any, onChanged: any = null) {
    this.dart = dart;
    if (onChanged != null)
      this.onChanged.subscribe((_) => onChanged(this.value));
  }

  /** Creates input for the specified property, and optionally binds it to the specified object */
  static forProperty(property: Property, source: any = null): InputBase {
    return toJs(api.grok_InputBase_ForProperty(property.dart, source));
  }

  /** Creates input for the specified input type */
  static forInputType(inputType: string | d4.InputType): InputBase {
    return toJs(api.grok_InputBase_ForInputType(inputType as string));
  }

  /** Creates input for the specified column */
  static forColumn<T = any>(column: Column<T>): InputBase<T | null> {
    return toJs(api.grok_InputBase_ForColumn(column.dart));
  }

  /** Input type identifier (such as "Slider" for the slider input). See {@link InputType}. */
  get inputType(): string {
    return api.grok_InputBase_Get_InputType(this.dart);
  };

  /** Data type this input can edit. See {@link Type}. */
  get dataType(): string { return api.grok_InputBase_Get_DataType(this.dart); };

  /** Visual root (typically a div element that contains {@link caption} and {@link input}) */
  get root(): HTMLElement { return api.grok_InputBase_Get_Root(this.dart); };

  get caption(): string {
    return api.grok_InputBase_Get_Caption(this.dart);
  }
  set caption(s: string) { api.grok_InputBase_Set_Caption(this.dart, s); }

  /** Property if associated with */
  get property(): any { return toJs(api.grok_InputBase_Get_Property(this.dart)); }
  set property(p: Property) { api.grok_InputBase_Set_Property(this.dart, toDart(p)); }

  /** Value format. */
  get format(): string { return api.grok_InputBase_Get_Format(this.dart); }
  set format(s: string) { api.grok_InputBase_Set_Format(this.dart, s); }

  get captionLabel(): HTMLElement { return api.grok_InputBase_Get_CaptionLabel(this.dart); }

  /** Returns the actual input */
  get input(): HTMLElement { return api.grok_InputBase_Get_Input(this.dart); }

  /** Whether empty values are allowed */
  get nullable(): boolean { return api.grok_InputBase_Get_Nullable(this.dart); }
  set nullable(v: boolean) { api.grok_InputBase_Set_Nullable(this.dart, v); }

  /** Whether events are thrown on value set */
  get notify(): boolean { return api.grok_InputBase_Get_Notify(this.dart); }
  set notify(v: boolean) { api.grok_InputBase_Set_Notify(this.dart, v); }

  /** Input value */
  get value(): T { return toJs(api.grok_InputBase_Get_Value(this.dart)); }
  set value(x: T) { api.grok_InputBase_Set_Value(this.dart, toDart(x)); }

  /** String representation of the {@link value} */
  get stringValue(): string { return api.grok_InputBase_Get_StringValue(this.dart); }
  set stringValue(s: string) { api.grok_InputBase_Set_StringValue(this.dart, s); }

  /** Whether the input is readonly */
  get readOnly(): boolean { return api.grok_InputBase_Get_ReadOnly(this.dart); }
  set readOnly(v: boolean) { api.grok_InputBase_Set_ReadOnly(this.dart, v); }

  /** Whether the input is enabled */
  get enabled(): boolean { return api.grok_InputBase_Get_Enabled(this.dart); }
  set enabled(v: boolean) { api.grok_InputBase_Set_Enabled(this.dart, v); }

  /** Occurs when [value] is changed, either by user or programmatically. */
  get onChanged(): Observable<T> { return api.grok_InputBase_OnChanged(this.dart); }

  /** Occurs when [value] is changed by user. */
  get onInput(): Observable<Event> { return api.grok_InputBase_OnInput(this.dart); }

  /** Saves the value. Used in dialog history. See also {@link load} */
  save(): any {
    return api.grok_InputBase_Save(this.dart);
  };

  /** Loads the value. Used in dialog history. See also {@link load} */
  load(s: any): any { return api.grok_InputBase_Load(this.dart, s); };

  init(): any {
    return api.grok_InputBase_Init(this.dart);
  };

  /** Fires the 'changed' event (value has changed). See also {@link fireInput} */
  fireChanged(): any {
    return api.grok_InputBase_FireChanged(this.dart);
  };

  /** Fires the 'input' event (user input). See also {@link fireChanged} */
  fireInput(): any {
    return api.grok_InputBase_FireInput(this.dart);
  };

  /** Adds the specified caption */
  addCaption(caption: string): InputBase<T> {
    api.grok_InputBase_AddCaption(this.dart, caption);
    return this;
  };

  /** Adds the specified postfix */
  addPostfix(postfix: string): InputBase<T> {
    api.grok_InputBase_AddPostfix(this.dart, postfix);
    return this;
  };

  /** Adds the specified options */
  addOptions(options: HTMLElement): InputBase<T> {
    api.grok_InputBase_AddOptions(this.dart, options);
    return this;
  };

  /** Adds a usage example to the input's hamburger menu */
  addPatternMenu(pattern: any): void {
    api.grok_InputBase_AddPatternMenu(this.dart, pattern);
  }

  /** Adds a validator that accepts a string representation of the edited value
   * and returns null if valid, or error message if invalid*/
  addValidator(validator: (value: string) => string | null): void {
    api.grok_InputBase_AddValidator(this.dart, validator);
  }

  /** Sets the tooltip */
  setTooltip(msg: string, tooltipCheck: (() => boolean) | null = null): InputBase<T> {
    api.grok_InputBase_SetTooltip(this.dart, msg, tooltipCheck);
    return this;
  };

  get classList(): DOMTokenList { return this.root.classList; }

  /**
   * Performs immediate validation of the input and returns the result.
   *
   * @returns {boolean} True if the input is valid; otherwise, false.
   */
  validate(): boolean {
    return api.grok_InputBase_Validate(this.dart);
  }
}


/** Base class for JS value editors */
export abstract class JsInputBase<T = any> extends InputBase<T> {
  //onInput: rxjs.Subject<any> = new rxjs.Subject<any>();

  abstract get inputType(): string;
  abstract get dataType(): string;

  abstract getInput(): HTMLElement;

  abstract getValue(): T;
  abstract setValue(value: T): void;

  abstract getStringValue(): string;
  abstract setStringValue(value: string): void;

  get input() { return this.getInput(); }

  get value(): T { return this.getValue(); }
  set value(value: T) { this.setValue(value); }

  get stringValue(): string { return this.getStringValue(); }
  set stringValue(value: string) { this.setStringValue(value); }

  constructor() {
    super(null, null);
    this.dart = api.grok_InputBase_FromJS(this);
  }
}

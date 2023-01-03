import {toDart, toJs} from "./wrappers";
import {__obs, _sub, observeStream, StreamSubscription} from "./events";
import {Observable, Subscription} from "rxjs";
import {Func, Property, PropertyOptions} from "./entities";
import {Cell, Column, DataFrame} from "./dataframe";
import {ColorType, FILTER_TYPE, LegendPosition, Type} from "./const";
import * as rxjs from "rxjs";
import { filter } from 'rxjs/operators';
import $ from "cash-dom";
import {MapProxy} from "./utils";
import dayjs from "dayjs";

declare let grok: any;
declare let DG: any;
declare let ui: any;
let api = <any>window;

export type RangeSliderStyle = 'barbell' | 'lines' | 'thin_barbell';

export type SliderOptions = {
  style?: RangeSliderStyle
}

export class ObjectPropertyBag {
  source: any;

  constructor(source: any, x: any = null) {

    /** @member {Object} */
    this.source = source;

    if (x == null)
      x = source;

    let props = this;

    const handler = {
      ownKeys(target: any) {
        return props.getProperties().map((p: Property) => p.name);
      },
      has(target: any, name: string) {
        return props.getProperties().find((p: Property) => p.name === name) !== null;
      },
      getOwnPropertyDescriptor(target: any, key: any) {
        return {
          enumerable: true,
          configurable: true
        };
      },
      set(target: any, name: string, value: any) {
        if (<any>x instanceof DartWidget)
          x = toDart(x);
        target.getProperty(name).set(x, value);
        return true;
      },
      get(target: any, name: any) {
        const own = ['__proto__', 'hasProperty', 'getProperty', 'getProperties', 'get', 'set', 'setAll', 'apply'];
        if (own.includes(name) || props.hasOwnProperty(name))
          // @ts-ignore
          return props[name];
        if (<any>x instanceof DartWidget)
          x = toDart(x);
        return target.getProperty(name)?.get(x);
      }
    }

    return new Proxy(this, handler);
  }

  /**
   * Gets the value of the specified property
   * @param {string} propertyName
   * @returns {object}
   * */
  get(propertyName: string): object {
    return this.getProperty(propertyName).get(this.source);
  }

  /**
   * Sets the value of the specified property
   * @param {string} propertyName
   * @param {object} propertyValue
   * */
  set(propertyName: string, propertyValue: object) {
    this.getProperty(propertyName).set(this.source, propertyValue);
  }

  /** Sets all properties according to the passed object containing key-value pairs */
  setAll(params: object) {
    for (const [k, v] of Object.entries(params))
      this.set(k, v);
  }

  /** @returns {Property[]} */
  getProperties(): Property[] {
    return this.source.getProperties();
  }

  /** Gets property by name (case-sensitive).
   * @param {string} name
   * @returns {Property} */
  getProperty(name: string): Property {
    let property = this.getProperties().find((p) => p.name === name);
    if (typeof property == 'undefined')
      throw `Property not found: ${name}`;
    return property;
  }

  /**
   * @param {string} name
   * @returns {boolean} */
  hasProperty(name: string): boolean {
    return this.getProperties().findIndex((p) => p.name === name) !== -1;
  }
}


/** Base class for controls that have a visual root and a set of properties. */
export class Widget<TSettings = any> {

  public get type(): string { return 'Unknown'; }

  /** Contains auxiliary information */
  public temp: any;

  /** Constructor function. No parameters, returns [Widget]. */
  factory: Func | null = null;

  protected _root: HTMLElement;
  protected _properties: Property[];
  props: TSettings & ObjectPropertyBag; //ObjectPropertyBag;
  subs: Subscription[];
  dart: any;
  isDetached: boolean = false;

  /** @constructs Widget and initializes its root. */
  constructor(widgetRoot: HTMLElement) {

    /** @member {HTMLElement} */
    this._root = widgetRoot;

    /** @member {Property[]}*/
    this._properties = [];

    /** @member {ObjectPropertyBag} */
    // @ts-ignore
    this.props = new ObjectPropertyBag(this);

    /** @member {StreamSubscription[]} */
    this.subs = [];

    this.getProperties = this.getProperties.bind(this);

    this.temp = {};
  }

  /** Returns all currently active widgets. */
  static getAll(): Widget[] {
    return api.grok_Widget_GetAll();
  }

  /** Finds existing widget from its visual root. */
  static find(root: Element): Widget | null {
    return api.grok_Widget_FromRoot(root);
  }

  toDart() {
    if (this.dart == null)
      this.dart = api.grok_Widget_Wrap(this);
    return this.dart;
  }

  /** Registers a subscription to an external event.
   * @param {Subscription} subscription */
  sub(subscription: Subscription): void {
    this.subs.push(subscription);
  }

  /**
   * @param {Object} properties
   * @returns {Widget} */
  apply(properties: object): Widget {
    for (let name of Object.keys(properties))
      if (typeof name !== 'undefined')
        // @ts-ignore
        this.props[name] = properties[name];

     //this.props.apply(properties);
    return this;
  }

  getProperties(): Property[] { return this._properties; }

  /** Gets called when viewer's property is changed.
   * @param {Property} property - or null, if multiple properties were changed. */
  onPropertyChanged(property: Property | null): void {}

  getDartProperties(): any[] {
    return this.getProperties().map((p) => p.dart);
  }

  onFrameAttached(dataFrame: DataFrame): void {
    if (this.props.hasProperty('dataFrame'))
      this.props.set('dataFrame', dataFrame);
  }

  /** Widget's visual root.
   * @type {HTMLElement} */
  get root(): HTMLElement { return this._root; }
  set root(r: HTMLElement) { this._root = r; }

  /** Gets called when a widget is detached and will no longer be used. Typically used for unsubscribing from events.
   * Be sure to call super.detach() if this method is overridden.  */
  detach(): void {
    this.subs.forEach((s) => s.unsubscribe());
    this.isDetached = true;
  }

  /** Registers an property with the specified type, name, and defaultValue.
   *  Registered property gets added to {@see properties}.
   *  Returns default value, thus allowing to combine registering a property with the initialization
   *
   * @param {string} propertyName
   * @param {TYPE} propertyType
   * @param defaultValue
   * @param {Object} options
   * @returns {*}
   * @private
   */
  addProperty(propertyName: string, propertyType: Type, defaultValue: any = null, options: { [key: string]: any } & PropertyOptions | null = null): any {
    const fieldName = options?.fieldName ?? propertyName;

    let obj = this;
    // @ts-ignore
    let p = Property.create(propertyName, propertyType, (_) => obj[fieldName], null, defaultValue);
    p.set = function (_: any, x: any) {
      // @ts-ignore
      obj[fieldName] = x;
      obj.onPropertyChanged(p);
    };

    if (options !== null) {
      for (let key of Object.keys(options))
        if (key != 'fieldName')
          api.grok_PropMixin_SetPropertyValue(p.dart, key, options[key]);
    }

    this._properties.push(p);
    return p.defaultValue;
  }

  /** Creates a new widget from the root element. */
  static fromRoot(root: HTMLElement): Widget {
    return new Widget(root);
  }

  // /** Creates a {@see Widget} from the specified React component. */
  // // @ts-ignore
  // static react(reactComponent: React.DOMElement<any, any> | Array<React.DOMElement<any, any>> | React.CElement<any, any> | Array<React.CElement<any, any>> | React.ReactElement | Array<React.ReactElement>): Widget {
  //   let widget = Widget.fromRoot(ui.div());
  //   // @ts-ignore
  //   ReactDOM.render(reactComponent, widget.root);
  //   return widget;
  // }
}


/** Base class for DataFrame-bound filtering controls.
 * Supports collaborative filtering by efficiently working together with
 * other filters. */
export abstract class Filter extends Widget {

  /** An indicator icon on the left of the filter header, before the name.
   * A filter is responsible for hiding or showing it, depending on its state. */
  indicator: HTMLDivElement;

  /** Group of control icons on the right of the filter header, after the name.
   * FilterGroup takes care of its visibility. */
  controls: HTMLDivElement;

  /** A DataFrame this filter is associated with. */
  dataFrame: DataFrame | null;

  /** A column this filter is associated with. */
  column: Column | null = null;

  columnName?: string;

  /** Caption to be shown in the filter group. */
  get caption(): string { return this.columnName ?? ''}

  constructor() {
    super(ui.div());

    this.dataFrame = null;
    this.indicator = ui.div([], 'd4-filter-indicator');
    this.controls = ui.div([], 'd4-flex-row');

    $(this.indicator).hide();
  }

  // static create(type: FILTER_TYPE | null, column: Column): Filter {
  //   return api.grok_Filter_ForColumn(column, type);
  // }
  //
  // static forColumn(column: Column) { return Filter.create(null, column); }
  //
  // static histogram(column: Column) { return Filter.create(FILTER_TYPE.HISTOGRAM, column); }
  // static categorical(column: Column) { return Filter.create(FILTER_TYPE.CATEGORICAL, column); }
  // static multiValue(column: Column) { return Filter.create(FILTER_TYPE.MULTI_VALUE, column); }
  // static freeText(column: Column) { return Filter.create(FILTER_TYPE.FREE_TEXT, column); }

  /** Override to indicate whether the filter actually filters something (most don't in the initial state).
   * This is used to minimize the number of unnecessary computations.
   * Make sure to call super.isFiltering to check whether the filter has been disabled by the user */
  get isFiltering(): boolean {
    return !(this.root.parentElement?.classList?.contains('d4-filter-disabled') == true) &&
      !(this.root.parentElement?.parentElement?.parentElement?.classList?.contains('d4-filters-disabled') == true);
  }

  /** Whether a filter is ready to apply the filtering mask synchronously. */
  get isReadyToApplyFilter(): boolean { return true; }

  /** Override to provide short filter summary that might be shown on viewers or in the property panel. */
  abstract get filterSummary(): string;

  /** Override to filter the dataframe.
   * The method should work with `this.dataFrame.filter`, should disregard
   * false values (these are filtered out already by other filters), and should
   * filter out corresponding indexes. */
  abstract applyFilter(): void;

  /** Override to save filter state. */
  saveState(): any {
    console.log('save state');
    return {
      column: this.columnName,
      columnName: this.columnName
    };
  }

  /** Override to load filter state. */
  applyState(state: any): void {
    this.columnName = state.columnName;
    this.column = this.columnName && this.dataFrame ? this.dataFrame.col(this.columnName) : null;
    console.log('apply state');
  }

  /** Gets called when a data frame is attached.
   * Make sure to call super.attach(dataFrame) when overriding.
   * @param {DataFrame} dataFrame*/
  attach(dataFrame: DataFrame): void {
    this.dataFrame = dataFrame;
    this.subs.push(this.dataFrame.onRowsFiltering
      .pipe(filter((_) => this.isFiltering && this.isReadyToApplyFilter))
      .subscribe((_) => {
        this.applyFilter();
        this.dataFrame!.rows.filters.push(`${this.columnName}: ${this.filterSummary}`);
      })
    );
  }

  /** Gets called when a previously used filter gets moved in the DOM tree.
   * Normally, you don't need to do anything, but this might be handy for
   * the iframe-based filters. */
  refresh() {}

  detach(): void {
    super.detach();
    if (this.isFiltering)
      this.dataFrame?.rows?.requestFilter();
  }
}


export class DartWidget extends Widget {
  constructor(dart: any) {
    super(api.grok_Widget_Get_Root(dart));
    this.dart = dart;
    this.temp = new MapProxy(api.grok_Widget_Get_Temp(this.dart));
  }

  get type(): string {
    return api.grok_Widget_Get_Type(this.dart);
  }

  get root(): HTMLElement {
    return api.grok_Widget_Get_Root(this.dart);
  }

  getProperties(): Property[] {
    return toJs(api.grok_PropMixin_GetProperties(this.dart));
  }
}

/**
 * Accordion control with collapsible/expandable panes.
 * Samples: {@link https://public.datagrok.ai/js/samples/ui/accordion}
 * @extends {DartWidget}
 * */
export class Accordion extends DartWidget {

  /** @constructs Accordion */
  constructor(dart: any) {
    super(dart);
  }

  /** An object this accordion is associated with */
  get context(): any {
    return toJs(api.grok_Accordion_Get_Context(this.dart));
  }
  set context(x: any) {
    api.grok_Accordion_Set_Context(this.dart, toDart(x));
  }

  /** Creates a new instance of Accordion */
  static create(key: any = null): Accordion {
    return toJs(api.grok_Accordion(key));
  }

  /** @type {AccordionPane[]} */
  get panes(): AccordionPane[] {
    return api.grok_TabControlBase_Get_Panes(this.dart).map(toJs);
  }

  /** Header element on top of the accordion */
  get header(): HTMLElement { return api.grok_Accordion_Get_Header(this.dart); }
  set header(header) { api.grok_Accordion_Set_Header(this.dart, header); }

  /** Returns a pane with the specified name.
   * @param {string} name
   * @returns {AccordionPane} */
  getPane(name: string): AccordionPane {
    return toJs(api.grok_TabControlBase_GetPane(this.dart, name));
  }

  /** Adds a title element. */
  addTitle(element: HTMLElement): void {
    return api.grok_Accordion_AddTitle(this.dart, element);
  }

  /** Adds a pane */
  addPane(name: string, getContent: () => HTMLElement, expanded: boolean = false, before: AccordionPane | null = null): AccordionPane {
    return toJs(api.grok_Accordion_AddPane(this.dart, name, getContent, expanded, before !== null ? before.dart : null, null));
  }

  /** Adds a pane with the count indicator next to the title.
   * getCount() is executed immediately. */
  addCountPane(name: string, getContent: () => HTMLElement, getCount: () => number, expanded: boolean = false, before: AccordionPane | null = null): AccordionPane {
    return toJs(api.grok_Accordion_AddPane(this.dart, name, getContent, expanded, before !== null ? before.dart : null, getCount));
  }

  /** Removed the specified pane. */
  removePane(pane: AccordionPane) {
    api.grok_Accordion_RemovePane(this.dart, pane.dart);
  }

  /** Finalizes accordion construction */
  end() {
    api.grok_Accordion_End(this.dart);
  }
}


/** A pane in the {@link Accordion} control. */
export class AccordionPane extends DartWidget {
  dart: any;

  constructor(dart: any) {
    super(dart);
  }

  /** Expanded state
   * @type {boolean} */
  get expanded(): boolean {
    return api.grok_AccordionPane_Get_Expanded(this.dart);
  }

  set expanded(v: boolean) {
    api.grok_AccordionPane_Set_Expanded(this.dart, v);
  }

  /** @type {string} */
  get name(): string {
    return api.grok_AccordionPane_Get_Name(this.dart);
  }

  set name(name: string) {
    api.grok_AccordionPane_Set_Name(this.dart, name);
  }
}


/** Tab control that hosts panes inside. See also {@link TabPane} */
export class TabControl {
  dart: any;
  constructor(dart: any) {
    this.dart = dart;
  }

  /** Creates a new TabControl */
  static create(vertical: boolean = false): TabControl {
    return toJs(api.grok_TabControl(vertical));
  }

  /** Visual root */
  get root(): HTMLDivElement {
    return api.grok_Widget_Get_Root(this.dart);
  }

  /** Header shown on top of the control */
  get header(): HTMLDivElement {
    return api.grok_TabControlBase_Get_Header(this.dart);
  }

  /** Panes currently present in the pane control.
   * Do not change the array, use {@link addPane} instead */
  get panes(): TabPane[] {
    return api.grok_TabControlBase_Get_Panes(this.dart).map(toJs);
  }

  /** Gets the pane with the specified name */
  getPane(name: string): TabPane {
    return toJs(api.grok_TabControlBase_GetPane(this.dart, name));
  }

  /** Adds a new pane with the specified name */
  addPane(name: string, getContent: () => HTMLElement, icon: any = null, options?: {allowClose: boolean}): TabPane {
    return toJs(api.grok_TabControlBase_AddPane(this.dart, name, getContent, icon, options?.allowClose ?? false));
  }

  /** Removes all panes */
  clear(): void {
    api.grok_TabControlBase_Clear(this.dart);
  }

  /** Currently visible pane */
  get currentPane(): TabPane { return api.grok_TabControlBase_Get_CurrentPane(this.dart); }
  set currentPane(v: TabPane) { api.grok_TabControlBase_Set_CurrentPane(this.dart, v.dart); }

  /** Occurs before the active pane is changed */
  get onBeforeTabChanged(): Observable<any> { return __obs('d4-tabcontrol-before-tab-changed', this.dart); }

  /** Occurs after the active pane is changed */
  get onTabChanged(): Observable<any> { return __obs('d4-tabcontrol-tab-changed', this.dart); }

  get onTabAdded(): Observable<any> { return __obs('d4-tabcontrol-tab-added', this.dart); }
  get onTabRemoved(): Observable<any> { return __obs('d4-tabcontrol-tab-removed', this.dart); }
}


/** Represents a pane of either {@link TabControl} or {@link Accordion} */
export class TabPane {
  dart: any;

  /** Creates TabPane from the Dart handle */
  constructor(dart: any) {
    this.dart = dart;
  }

  /** {@link TabControl} this pane belongs to */
  get parent(): TabControl {
    return toJs(api.grok_TabPane_Get_Parent(this.dart));
  }

  /** A control shown on top of the pane */
  get header(): HTMLDivElement {
    return api.grok_TabPane_Get_Header(this.dart);
  }

  /** Content */
  get content(): HTMLDivElement {
    return api.grok_TabPane_Get_Content(this.dart);
  }

  /** Whether the pane is expanded. Applicable to Accordion's panes only. */
  get expanded(): boolean { return api.grok_AccordionPane_Get_Expanded(this.dart); }
  set expanded(v: boolean) { api.grok_AccordionPane_Set_Expanded(this.dart, v); }

  /** Tab pane name */
  get name(): string { return api.grok_AccordionPane_Get_Name(this.dart); }
  set name(name: string) { api.grok_AccordionPane_Set_Name(this.dart, name); }
}


export class ToolboxPage {
  dart: any;
  constructor(dart: any) {
    this.dart = dart;
  }

  get accordion(): Accordion {
    return toJs(api.grok_ToolboxPage_Get_Accordion(this.dart));
  }
}


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
export class Dialog extends DartWidget {

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
   * @returns {Dialog} */
  onOK(handler: Function): Dialog {
    api.grok_Dialog_OnOK(this.dart, handler);
    return this;
  }

  /**
   * Sets the CANCEL button handler
   * @param {Function} handler
   * @returns {Dialog} */
  onCancel(handler: Function): Dialog {
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
  show(options?: { modal?: boolean; resizable?: boolean; fullScreen?: boolean; center?: boolean; centerAt?: Element; x?: number; y?: number; width?: number; height?: number; backgroundColor?: string; }): Dialog {
    api.grok_Dialog_Show(this.dart, options?.modal, options?.resizable, options?.fullScreen, options?.center, options?.centerAt, options?.x, options?.y, options?.width, options?.height, options?.backgroundColor);
    return this;
  }

  /** @returns {Dialog}
   * @param {boolean} fullScreen  */
  showModal(fullScreen: boolean): Dialog {
    api.grok_Dialog_Show(this.dart, true, null, fullScreen, false, null, null, null, null, null, null);
    return this;
  }

  /** Adds content to the dialog.
   * @param {HTMLElement | Widget | InputBase} content
   * @returns {Dialog} */
  add(content: HTMLElement | Widget | InputBase): Dialog {
    api.grok_Dialog_Add(this.dart, toDart(content));
    return this;
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
  addButton(text: string, action: Function, index: number = 0, tooltip: any = null): Dialog {
    api.grok_Dialog_AddButton(this.dart, text, action, index, tooltip);
    return this;
  }

  /** Adds context action with the specified text.
   * @param {string} text
   * @param {Function} action
   * @returns {Dialog}
   * */
  addContextAction(text: string, action: Function): Dialog {
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
  initDefaultHistory(): Dialog {
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
}


/** See {@link Menu.items} */
export interface IMenuItemsOptions<T = any> {

  /** Whether or not a check box appears before the item */
  isChecked?: (item: T) => boolean;

  /** If result is not null, the item is grayed out and the result is shown in the tooltip */
  isValid?: (item: T) => string | null;

  /** Text to be shown on the menu item */
  toString?: (item: T) => string;

  /** Tooltip */
  getTooltip?: (item: T) => string;

  /** Gets invoked when the mouse enters the item */
  onMouseEnter?: (item: T) => void;
}


/**
 * Menu (either top menu or popup menu).
 * Top menu sample: {@link https://public.datagrok.ai/js/samples/ui/menu}
 * Popup menu sample: {@link https://public.datagrok.ai/js/samples/ui/popup-menu}
 *
 * @example
 * DG.Menu.popup()
 *   .item('Show info', () => grok.shell.info('Info'))
 *   .separator()
 *   .items(['First', 'Second'], showBalloon)
 *   .show();
 * */
export class Menu {
  dart: any;
  _check: HTMLDivElement = ui.div('d4-menu-item-check');

  constructor(dart: any) {
    this.dart = dart;
  }

  static create(): Menu {
    return toJs(api.grok_Menu());
  }

  /** Creates a popup menu. */
  static popup(): Menu {
    return toJs(api.grok_Menu_Context());
  }

  get root(): HTMLElement { return api.grok_Menu_Get_Root(this.dart); }

  /** Finds a child menu item with the specified text. */
  find(text: string): Menu {
    return toJs(api.grok_Menu_Find(this.dart, text));
  }

  /** Executes the onClick function for that menu item.
   * Only works for items, not groups. */
  click(): void {
    api.grok_Menu_Click(this.dart);
  }

  /** Removes a child menu item with the specified text. */
  remove(text: string): void {
    api.grok_Menu_Remove(this.dart, text);
  }

  /** Removes all child menu items. */
  clear(): void {
    api.grok_Menu_Clear(this.dart);
  }

  /** Returns an existing menu group or adds a new group with the specified text. */
  group(text: string, order: number | null = null): Menu {
    return toJs(api.grok_Menu_Group(this.dart, text, order));
  }

  /** Ends a group of menu items and returns to the higher menu level.
   * @returns {Menu} */
  endGroup(): Menu {
    return toJs(api.grok_Menu_EndGroup(this.dart));
  }

  /** Adds a menu group with the specified text and handler. */
  item(text: string, onClick: () => void, order: number | null = null): Menu {
    return toJs(api.grok_Menu_Item(this.dart, text, onClick, order));
  }

  /** For each item in items, adds a menu group with the specified text and handler. */
  items<T = any>(items: T[], onClick: (item: T) => void, options: IMenuItemsOptions<T> | null = null): Menu {
    return toJs(api.grok_Menu_Items(this.dart, items, onClick, options?.isValid, options?.isChecked, options?.toString, options?.getTooltip, options?.onMouseEnter));
  }

  /** Adds a separator line.
   *  @returns {Menu} */
  separator(): Menu {
    return toJs(api.grok_Menu_Separator(this.dart));
  }

  /** Shows the menu.
   * @returns {Menu} */
  show(): Menu {
    return toJs(api.grok_Menu_Show(this.dart));
  }

  get onContextMenuItemClick() {
    return __obs('d4-menu-item-click', this.dart);
  }

  toString(): string {
    return api.grok_MenuItem_ToString(this.dart);
  }
}


/** Balloon-style visual notifications. */
export class Balloon {

  /** Shows information message (green background) */
  info(s: string): void {
    api.grok_Balloon(s, 'info');
  }

  /** Shows information message (red background) */
  error(s: string): void {
    api.grok_Balloon(s, 'error');
  }

  /** Closes all balloons currently shown */
  static closeAll(): void {
    api.grok_Balloon_CloseAll();
  }
}



/** Input control base. Could be used for editing {@link Property} values as well.
 * The root is a div that consists of {@link captionLabel} and {@link input}.
 * */
export class InputBase<T = any> {
  dart: any;

  constructor(dart: any, onChanged: any = null) {
    this.dart = dart;
    if (onChanged != null)
      this.onChanged((_: any) => onChanged(this.value));
  }

  /** Creates input for the specified property, and optionally binds it to the specified object */
  static forProperty(property: Property, source: any = null): InputBase {
    return toJs(api.grok_InputBase_ForProperty(property.dart, source));
  }

  /** Creates input for the specified column */
  static forColumn<T = any>(column: Column<T>): InputBase<T | null> {
    return toJs(api.grok_InputBase_ForColumn(column.dart));
  }

  /** Visual root (typically a div element that contains {@link caption} and {@link input}) */
  get root(): HTMLElement { return api.grok_InputBase_Get_Root(this.dart); };

  get caption(): string {
    return api.grok_InputBase_Get_Caption(this.dart);
  }

  /** Property if associated with */
  get property(): any { return toJs(api.grok_InputBase_Get_Property(this.dart)); }

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
  set value(x: T) { toDart(api.grok_InputBase_Set_Value(this.dart, x)); }

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
  onChanged(callback: Function): StreamSubscription {
    return _sub(api.grok_InputBase_OnChanged(this.dart, callback));
  }

  /** Occurs when [value] is changed by user. */
  onInput(callback: Function): StreamSubscription {
    return _sub(api.grok_InputBase_OnInput(this.dart, callback));
  }

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

  /** Sets the tooltip */
  setTooltip(msg: string, tooltipCheck: (() => boolean) | null = null): InputBase<T> {
    api.grok_InputBase_SetTooltip(this.dart, msg, tooltipCheck);
    return this;
  };

  get classList(): DOMTokenList { return this.root.classList; }
}


/** Base class for JS value editors */
export abstract class JsInputBase<T = any> extends InputBase<T> {
  //onInput: rxjs.Subject<any> = new rxjs.Subject<any>();

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
    //this.dart = api.grok_InputBase_FromJS(this);
  }
}


export class DateInput extends InputBase<dayjs.Dayjs> {
  dart: any;

  constructor(dart: any, onChanged: any = null) {
    super(dart, onChanged);
  }

  get value(): dayjs.Dayjs { return dayjs(api.grok_DateInput_Get_Value(this.dart)); }
  set value(x: dayjs.Dayjs) { toDart(api.grok_DateInput_Set_Value(this.dart, x.valueOf())); }
}

export class ProgressIndicator {
  dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  static create() {
    return toJs(api.grok_ProgressIndicator_Create());
  }

  get percent(): number {
    return api.grok_ProgressIndicator_Get_Percent(this.dart);
  }

  get description(): string {
    return api.grok_ProgressIndicator_Get_Description(this.dart);
  }

  set description(s: string) {
    api.grok_ProgressIndicator_Set_Description(this.dart, s);
  }

  update(percent: number, description: string): void {
    api.grok_ProgressIndicator_Update(this.dart, percent, description);
  }

  log(line: string): void {
    api.grok_ProgressIndicator_Log(this.dart, line);
  }

  get onProgressUpdated(): Observable<any> {
    return observeStream(api.grok_Progress_Updated(this.dart));
  }

  get onLogUpdated(): Observable<any> {
    return observeStream(api.grok_Progress_Log_Updated(this.dart));
  }
}


export class TaskBarProgressIndicator extends ProgressIndicator {
  static create(name?: string): TaskBarProgressIndicator {
    return toJs(api.grok_TaskBarProgressIndicator_Create(name));
  }

  close(): any {
    return api.grok_TaskBarProgressIndicator_Close(this.dart);
  }
}


export class TagEditor {
  dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  static create(): TagEditor {
    return toJs(api.grok_TagEditor());
  }

  get root(): HTMLElement {
    return api.grok_TagEditor_Get_Root(this.dart);
  }

  get tags(): TagElement[] {
    return api.grok_TagEditor_Get_Tags(this.dart);
  }

  addTag(tag: TagElement | string, notify: boolean = true) {
    return api.grok_TagEditor_AddTag(this.dart, tag, notify);
  }

  removeTag(tag: TagElement | string): void {
    api.grok_TagEditor_RemoveTag(this.dart, tag);
  }

  clearTags(): void {
    api.grok_TagEditor_ClearTags(this.dart);
  }

  set acceptsDragDrop(predicate: (...params: any[]) => boolean) {
    api.grok_TagEditor_Set_AcceptsDragDrop(this.dart, (x: any) => predicate(toJs(x, false)));
  };

  set doDrop(action: Function) {
    api.grok_TagEditor_Set_DoDrop(this.dart, (x: any) => action(toJs(x, false)));
  }

  onChanged(callback: Function): StreamSubscription {
    return _sub(api.grok_TagEditor_OnChanged(this.dart, callback));
  }
}


export class TagElement {
  dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  get tag(): TagElement | string {
    return api.grok_TagElement_Get_Tag(this.dart);
  };

  set tag(x: TagElement | string) {
    api.grok_TagElement_Set_Tag(this.dart, x);
  };
}


/** Color-related routines. */
export class Color {

  /** Returns a color associated with the specified cell as an ARGB-formatted integer.
   * To convert to html color, use {@link getCellColorHtml} or {@link toHtml}. */
  static getCellColor(cell: Cell): number {
    return api.grok_Color_FromCell(cell.dart);
  }

  /** Returns a string representation of the color associated with the specified cell.
   * For batch color manipulations, use {@link getCellColor}. */
  static getCellColorHtml(cell: Cell): string {
    return Color.toHtml(Color.getCellColor(cell));
  }

  /** Returns a color associated with the specified category within a column.
   * Returns ARGB-formatted integer. To convert to html color, use {@link toHtml}. */
  static getCategoryColor(column: Column, category: any): number {
    return api.grok_Color_FromCategory(column.dart, category);
  }

  /** Returns the Alpha component of the color represented as ARGB-formatted integer. */
  static a(c: number): number { return (c >> 24) & 0xFF; }

  /** Returns the Red component of the color represented as ARGB-formatted integer. */
  static r(c: number): number { return (c >> 16) & 0xFF; }

  /** Returns the Green component of the color represented as ARGB-formatted integer. */
  static g(c: number): number { return (c >> 8) & 0xFF; }

  /** Returns the Blue component of the color represented as ARGB-formatted integer. */
  static b(c: number): number { return c & 0xFF; }

  /** Returns i-th categorical color (looping over the palette if needed) */
  static getCategoricalColor(i: number): ColorType {
    return Color.categoricalPalette[i % Color.categoricalPalette.length];
  }

  /** Returns either black or white color, depending on which one would be most contrast to the specified [color]. */
  static getContrastColor(color: ColorType): ColorType {
    return api.grok_Color_GetContrastColor(color);
  }

 /** Converts ARGB-formatted integer color to a HTML-formatted string (such as `#ffffff`). See also {@link toRgb}. */
  static toHtml(color: number): string { return api.grok_Color_ToHtml(color); }

  /** Convert HTML-formatted string (such as `#ffffff`) to ARGB-fromatted integer color */
  static fromHtml(htmlColor: string): number { return api.grok_Color_FromHtml(htmlColor); }

  /** Converts ARGB-formatted integer color to a HTML-formatted string (such as `rbg(20, 46, 124)`). See also {@link toHtml. }*/
  static toRgb(color: ColorType): string {
    return color === null ? '' : `rgb(${Color.r(color)},${Color.g(color)},${Color.b(color)})`;
  }

  /** Returns the standard palette of the categorical colors used across all visualizations in Datagrok. */
  static get categoricalPalette(): ColorType[] {
    return api.grok_Color_CategoricalPalette();
  }

  static scaleColor(x: number, min: number, max: number, alpha?: number, colorScheme?: number[]): number {
    return api.grok_Color_ScaleColor(x, min, max, alpha ? alpha : null, colorScheme ? colorScheme : null);
  }

  static scale(x: number, min: number, max: number): number {
    return min === max ? min : (x - min) / (max - min);
  }

  static get gray(): ColorType {
    return 0xFF808080;
  }

  static get lightLightGray(): ColorType {
    return 0xFFF0F0F0;
  }

  static get lightGray(): ColorType {
    return 0xFFD3D3D3;
  }

  static get darkGray(): ColorType {
    return 0xFF838383;
  }

  static get blue(): ColorType {
    return 0xFF0000FF;
  }

  static get green(): ColorType {
    return 0xFF00FF00;
  }

  static get darkGreen(): ColorType {
    return 0xFF006400;
  }

  static get black(): ColorType {
    return 0xFF000000;
  }

  static get yellow(): ColorType {
    return 0xFFFFFF00;
  }

  static get white(): ColorType {
    return 0xFFFFFFFF;
  }

  static get red(): ColorType {
    return 0xFFFF0000;
  }

  static get darkRed(): ColorType {
    return 0xFF8b0000;
  }

  static get maroon(): ColorType {
    return 0xFF800000;
  }

  static get olive(): ColorType {
    return 0xFF808000;
  }

  static get orange(): ColorType {
    return 0xFFFFA500;
  }

  static get darkOrange(): ColorType {
    return 0xFFFF8C00;
  }

  static get lightBlue(): ColorType {
    return 0xFFADD8E6;
  }

  static get darkBlue(): ColorType {
    return 0xFF0000A0;
  }

  static get purple(): ColorType {
    return 0xFF800080;
  }

  static get whitesmoke(): ColorType {
    return 0xFFF5F5F5;
  }

  static get navy(): ColorType {
    return 0xFF000080;
  }

  static get cyan(): ColorType {
    return 0xFF00ffff;
  }

  static get filteredRows(): ColorType {
    return 0xff1f77b4;
  }

  static get filteredOutRows(): ColorType {
    return Color.lightLightGray;
  }

  static get selectedRows(): ColorType {
    return Color.darkOrange;
  }

  static get missingValueRows(): ColorType {
    return Color.filteredOutRows;
  }

  static get mouseOverRows(): ColorType {
    return 0xFFAAAAAA;
  }

  static get currentRow(): ColorType {
    return 0xFF38B738;
  }

  static get histogramBar(): ColorType {
    return Color.filteredRows;
  }

  static get barChart(): ColorType {
    return 0xFF24A221;
  }

  static get scatterPlotMarker(): ColorType {
    return 0xFF40699c;
  }

  static get scatterPlotSelection(): ColorType {
    return 0x80323232;
  }

  static get scatterPlotZoom(): ColorType {
    return 0x80626200;
  }

  static get areaSelection(): ColorType {
    return Color.lightBlue;
  }

  static get rowSelection(): ColorType {
    return 0x60dcdca0;
  }

  static get colSelection(): ColorType {
    return 0x60dcdca0;
  }

  static get areaZoom(): ColorType {
    return 0x80323232;
  }

  static get gridWarningBackground(): ColorType {
    return 0xFFFFB9A7;
  }

  static get success(): ColorType {
    return 0xFF3cb173;
  }

  static get failure(): ColorType {
    return 0xFFeb6767;
  }
}


/** Tree view node.
 * Sample: {@link https://public.datagrok.ai/js/samples/ui/tree-view}
 * */
export class TreeViewNode {
  dart: any;

  /** @constructs {TreeView} from the Dart object */
  constructor(dart: any) {
    this.dart = dart;
  }

  /** Visual root */
  get root(): HTMLElement {
    return api.grok_TreeViewNode_Root(this.dart);
  }

  /* Node's parent */
  get parent(): TreeViewNode {
    return api.grok_TreeViewNode_Parent(this.dart);
  }

  /** Caption label */
  get captionLabel(): HTMLElement {
    return api.grok_TreeViewNode_CaptionLabel(this.dart);
  }

  /** Check box element  */
  get checkBox(): HTMLElement | null {
    return api.grok_TreeViewNode_CheckBox(this.dart);
  }

  /** Returns `true` if checked */
  get checked(): boolean { return api.grok_TreeViewNode_Get_Checked(this.dart); }
  set checked(checked: boolean) { api.grok_TreeViewNode_Set_Checked(this.dart, checked); }

  /** Node text */
  get text(): string { return api.grok_TreeViewNode_Text(this.dart); }

  /** Node value */
  get value(): object { return api.grok_TreeViewNode_Get_Value(this.dart); };
  set value(v: object) { api.grok_TreeViewNode_Set_Value(this.dart, v); };

  /** Enables checkbox */
  enableCheckBox(checked: boolean = false): void {
    api.grok_TreeViewNode_EnableCheckBox(this.dart, checked);
  }

  /**  */
  get onSelected(): Observable<TreeViewNode> { return __obs('d4-tree-view-node-current', this.dart); }

  /** Removes the node and its children from the parent */
  remove(): void {
    api.grok_TreeViewNode_Remove(this.dart);
  }
}

export class TreeViewGroup extends TreeViewNode {
  /** Creates new tree */
  static tree(): TreeViewGroup {
    return toJs(api.grok_TreeViewNode_Tree());
  }

  static fromItemCategories(items: any[], props: string[], options?: {
    removeEmpty: boolean, itemToElement?: (item:any) => Element, itemToString?: (item: any) => string, itemToValue?: (item: any) => any
  }): TreeViewGroup {

    function init(node: TreeViewGroup, path: string[]) {

      //
      const pathItems = items
        .filter((item) => props.every((p: string) => !(options?.removeEmpty ?? false) || item[p] != null))
        .filter((item) => path.every((p, i) => item[props[i]] == path[i]));

      // leafs
      if (path.length == props.length) {
        let itemToValue = options!.itemToValue ?? function (i) {return i;};
        let itemToString = options!.itemToElement ?? options!.itemToString;
        for (let item of pathItems)
          node.item(itemToString!(item), itemToValue(item));
      }
      else {
        let categories: Set<string> = new Set<string>();
        for (let item of pathItems)
          categories.add(item[props[path.length]]);

        for (let category of categories)
          init(node.group(category), [...path, category]);

        // lazy loading - use it for big collections and async queries
        // for (let category of categories)
        //   node.group(category)
        //     .onNodeExpanding
        //     .subscribe((expandingNode) => init(expandingNode, [...path, category]));
      }
    }

    let rootNode = TreeViewGroup.tree();
    init(rootNode, []);
    return rootNode;
  }

  /** Gets all node items */
  get items(): TreeViewNode[] {
    return api.grok_TreeViewNode_Items(this.dart).map((i: any) => toJs(i));
  }

  /** Gets the node's children */
  get children(): TreeViewNode[] {
    return api.grok_TreeViewNode_Children(this.dart).map((i: any) => toJs(i));
  }

  get expanded(): boolean { return api.grok_TreeViewNode_Get_Expanded(this.dart); }

  set expanded(isExpanded: boolean) { api.grok_TreeViewNode_Set_Expanded(this.dart, isExpanded); }

  /** Indicates whether check or uncheck is applied to a node only or to all node's children */
  get autoCheckChildren(): boolean { return api.grok_TreeViewNode_GetAutoCheckChildren(this.dart); }
  set autoCheckChildren(auto: boolean) { api.grok_TreeViewNode_SetAutoCheckChildren(this.dart, auto); }

  /** Adds new group */
  group(text: string | Element, value: object | null = null, expanded: boolean = true): TreeViewGroup {
    return toJs(api.grok_TreeViewNode_Group(this.dart, text, value, expanded));
  }

  /** Returns existing, or creates a new node group */
  getOrCreateGroup(text: string, value: object | null = null, expanded: boolean = true): TreeViewGroup {
    return toJs(api.grok_TreeViewNode_GetOrCreateGroup(this.dart, text, value, expanded));
  }

  /** Adds new item to group */
  item(text: string | Element, value: object | null = null): TreeViewNode {
    return toJs(api.grok_TreeViewNode_Item(this.dart, text, value));
  }

  get onNodeExpanding(): Observable<TreeViewGroup> { return __obs('d4-tree-view-node-expanding', this.dart); }

  get onNodeAdded(): Observable<TreeViewNode> { return __obs('d4-tree-view-node-added', this.dart); }

  get onNodeCheckBoxToggled(): Observable<TreeViewNode> { return __obs('d4-tree-view-node-checkbox-toggled', this.dart); }

  get onChildNodeExpandedChanged(): Observable<TreeViewGroup> { return __obs('d4-tree-view-child-node-expanded-changed', this.dart); }

  get onChildNodeExpanding(): Observable<TreeViewGroup> { return __obs('d4-tree-view-child-node-expanding', this.dart); }

  // get onChildNodeContextMenu(): Observable<TreeViewNode> { return __obs('d4-tree-view-child-node-context-menu', this.dart); }
  get onNodeContextMenu(): Observable<TreeViewNode> { return __obs('d4-tree-view-node-context-menu', this.dart); }

  get onSelectedNodeChanged(): Observable<TreeViewNode> { return __obs('d4-tree-view-selected-node-changed', this.dart); }

  get onNodeMouseEnter(): Observable<TreeViewNode> { return __obs('d4-tree-view-child-node-mouse-enter', this.dart); }

  get onNodeMouseLeave(): Observable<TreeViewNode> { return __obs('d4-tree-view-child-node-mouse-leave', this.dart); }

  get onNodeEnter(): Observable<TreeViewNode> { return __obs('d4-tree-view-node-enter', this.dart); }
}


/** A slider that lets user control both min and max values. */
export class RangeSlider extends DartWidget {

  static create(vertical: boolean = false, optionsStyle: RangeSliderStyle | SliderOptions = 'barbell'): RangeSlider {
    let options: SliderOptions = {};
    if (typeof optionsStyle !== 'string') {
      options.style = optionsStyle.style;
    } else {
      options.style = optionsStyle;
    }
    return toJs(api.grok_RangeSlider(options.style, vertical));
  }

  /** Minimum range value. */
  get minRange(): number { return api.grok_RangeSlider_Get_MinRange(this.dart); };

  /** Gets maximum range value. */
  get maxRange(): number { return api.grok_RangeSlider_Get_MaxRange(this.dart); };

  /** Gets minimum value. */
  get min(): number { return api.grok_RangeSlider_Get_Min(this.dart); };

  /** Gets maximum value. */
  get max(): number { return api.grok_RangeSlider_Get_Max(this.dart); };

  /** Sets values to range slider.
   * @param {number} minRange
   * @param {number} maxRange
   * @param {number} min
   * @param {number} max */
  setValues(minRange: number, maxRange: number, min: number, max: number): void {
    api.grok_RangeSlider_SetValues(this.dart, minRange, maxRange, min, max);
  }

  /** Sets showHandles in range slider.
   * @param {boolean} value */
  setShowHandles(value: boolean): void {
    console.log('setShowHandles', value);
    api.grok_RangeSlider_SetShowHandles(this.dart, value);
  }

  /** Sets the specified min value, preserving the range (scrolls). */
  scrollTo(newMinValue: number): void { return api.grok_RangeSlider_ScrollTo(this.dart, newMinValue); }

  /** Shifts min and max values by the specified delta. */
  scrollBy(delta: number): void { return api.grok_RangeSlider_ScrollTo(this.dart, delta); }

  /** @returns {Observable} */
  get onValuesChanged(): Observable<any> {
    return observeStream(api.grok_RangeSlider_Get_OnValuesChanged(this.dart));
  }
}


export class HtmlTable extends DartWidget {

  /** @constructs {HtmlTable} */
  constructor(dart: any) {
    super(dart);
  }

  /** Creates a visual table based on [items], [renderer], and [columnNames]
   * @returns {HtmlTable} */
  static create(items: any, renderer: any, columnNames: any = null): HtmlTable {
    return toJs(api.grok_HtmlTable(items, renderer !== null ? (object: any, ind: number) => renderer(toJs(object), ind) : null, columnNames));
  }

  /** Removes item */
  remove(item: any): void {
    api.grok_HtmlTable_Remove(this.dart, item);
  }
}


/** A combo box with columns as item.
 * Supports sorting, searching, custom tooltips, column-specific rendering, drag-and-drop, etc */
export class ColumnComboBox extends DartWidget {
  /** @constructs {ColumnComboBox} */
  constructor(dart: any) {
    super(dart);
  }

  /** Creates a column combo box with specified [dataframe] and [predicate]. */
  static create(dataframe: DataFrame, predicate: Function): ColumnComboBox {
    return toJs(api.grok_ColumnComboBox(dataframe.dart, (x: any) => predicate(toJs(x))));
  }

  /** @type {boolean} */
  get vertical(): boolean {
    return api.grok_ColumnComboBox_Get_Vertical(this.dart);
  }

  /** Text to be shown before the combo box */
  get caption(): string { return api.grok_ColumnComboBox_Get_Caption(this.dart); }
  set caption(c: string) { api.grok_ColumnComboBox_Set_Caption(this.dart, c); }

  /** @type {Property} */
  get property(): Property {
    return toJs(api.grok_ColumnComboBox_Get_Property(this.dart));
  }

  onEvent(eventId: string): rxjs.Observable<any> {
    return __obs(eventId, this.dart);
  }

  /** Occurs when the value is changed. */
  get onChanged(): rxjs.Observable<String> { return this.onEvent('d4-column-box-column-changed'); }
}


/** Column legend for viewers */
export class Legend extends DartWidget {
  constructor(dart: any) {
    super(dart);
  }

  static create(column: Column): Legend {
    return toJs(api.grok_Legend(column.dart));
  }

  /** Column for the legend */
  get column(): Column { return toJs(api.grok_Legend_Get_Column(this.dart)); }
  set column(column: Column) { api.grok_Legend_Set_Column(this.dart, column.dart); }

  /** Whether or not to show empty categories */
  get showNulls(): Boolean { return api.grok_Legend_Get_ShowNulls(this.dart); }
  set showNulls(show: Boolean) { api.grok_Legend_Set_ShowNulls(this.dart, show); }

  /** Position (left / right / top / bottom) */
  get position(): LegendPosition { return api.grok_Legend_Get_Position(this.dart); }
  set position(pos: LegendPosition) { api.grok_Legend_Set_Position(this.dart, pos); }
}


export class PropertyGrid extends DartWidget {

  constructor() {
    super(api.grok_PropertyGrid());
  }

  update(src: any, props: Property[]) {
    api.grok_PropertyGrid_Update(this.dart, src, props.map((x) => toDart(x)));
  }
}

export type fileShares = 'S3';

/** File browser widget */
export class FilesWidget extends DartWidget {
  /** Creates a [FilesWidget] and opens a directory, if [path] is specified.
   * [path] accepts a full-qualified name (see [Entity.nqName]). */
  static create(params: {path?: string, dataSourceFilter?: fileShares[]} = {}): FilesWidget {
    return toJs(api.grok_FilesWidget(params));
  }
}

/**
 * Routines for building UI
 * @module ui
 **/

import {ElementOptions, IndexPredicate} from './src/const';
import {Viewer} from './src/viewer';
import {VirtualView} from './src/view';
import {Accordion, Dialog, InputBase, Menu, TabControl, TreeViewNode, Widget, RangeSlider} from './src/widgets';
import {toDart, toJs} from './src/wrappers';
import {Functions} from './src/functions';
import $ from 'cash-dom';
import {__obs} from './src/events';
import {_isDartium, _options} from './src/utils';
import * as rxjs from 'rxjs';
import { CanvasRenderer, GridCellRenderer } from './src/grid';
import { DateTime } from './src/entities';
import { Column, DataFrame } from './src/dataframe';


let api = <any>window;

/**
 * Creates an instance of the element for the specified tag, and optionally assigns it a CSS class.
 * @param {string} tagName The name of an element.
 * @param {string | null} className
 * @returns {HTMLElement}
 */
export function element(tagName: string, className: string | null = null): HTMLElement {
  let x = document.createElement(tagName);
  if (className !== null)
    _class(x, className);
  return x;
}

/** Appends multiple elements to root, and returns root.
 *  An element could be either {@link HTMLElement} or {@link Viewer}.
 *
 * @param {HTMLElement} root
 * @param {(HTMLElement | Viewer)[]} elements
 * @returns {HTMLElement} */
export function appendAll(root: HTMLElement, elements: (HTMLElement | Viewer)[]): HTMLElement {
  let fragment = document.createDocumentFragment();
  for (let e of elements)
    fragment.appendChild(render(e));
  root.appendChild(fragment);
  return root;
}

export function empty(e: HTMLElement): HTMLElement {
  while (e.firstChild)
    e.removeChild(e.firstChild);
  return e;
}

export function _innerText(x: HTMLElement, s: string): HTMLElement {
  x.innerText = s;
  return x;
}

export function _class(x: HTMLElement, s: string): HTMLElement {
  x.classList.add(s);
  return x;
}

export function _color(x: HTMLElement, s: string): HTMLElement {
  x.style.color = s;
  return x;
}

export function _backColor(x: HTMLElement, s: string): HTMLElement {
  x.style.backgroundColor = s;
  return x;
}

/**
 * @param {number} height
 * @param {number} width
 * @returns {HTMLCanvasElement}
 * */
export function canvas(width: number | null = null, height: number | null = null): HTMLCanvasElement {
  let result = element('CANVAS');
  if (height != null && width != null) {
    $(result).height(height);
    $(result).width(width);
    $(result).prop('height', `${height}`);
    $(result).prop('width', `${width}`);
  }
  return result as HTMLCanvasElement;
}

/** @returns {HTMLHeadingElement} */
export function h1(s: string): HTMLHeadingElement {
  return _innerText(element('h1'), s) as HTMLHeadingElement;
}

/** @returns {HTMLHeadingElement} */
export function h2(s: string): HTMLHeadingElement {
  let x = element('h2');
  x.innerText = s;
  return x as HTMLHeadingElement;
}

/** @returns {HTMLHeadingElement} */
export function h3(s: string): HTMLHeadingElement {
  let x = element('h3');
  x.innerText = s;
  return x as HTMLHeadingElement;
}

/** @returns {Accordion} */
export function accordion(key: any = null): Accordion {
  return Accordion.create(key);
}

/**
 * @param {Object} pages - list of page factories
 * @param {boolean} vertical
 * @returns {TabControl} */
export function tabControl(pages: { [key: string]: any; } | null = null, vertical: boolean = false): TabControl {
  let tabs = TabControl.create(vertical);
  if (pages != null) {
    for (let key of Object.keys(pages)) {
      let value = pages[key];
      tabs.addPane(key, value instanceof Function ? value : () => render(value));
    }
  }
  return tabs;
}

/** Returns DivElement with the specified inner text
 * @param {string} text
 * @param {string | ElementOptions | null} options
 * @returns {HTMLDivElement} */
export function divText(text: string, options: string | ElementOptions | null = null): HTMLDivElement {
  let e = element('div');
  e.innerText = text;
  _options(e, options);
  return e as HTMLDivElement;
}

/** Returns a font-awesome icon with the specified name, handler, and tooltip.
 * Sample: {@link https://public.datagrok.ai/js/samples/ui/icons}
 * @param {string} name - icon name (omit the "fa-" prefix)
 * @param {Function} handler
 * @param {String} tooltipMsg
 * @returns {HTMLElement} */
export function iconFA(name: string, handler: (this: HTMLElement, ev: MouseEvent) => any, tooltipMsg: string | null = null): HTMLElement {
  let i = element('i');
  i.classList.add('grok-icon');
  i.classList.add('fal');
  i.classList.add(`fa-${name}`);
  i.addEventListener('click', handler);
  if (tooltipMsg !== null)
    tooltip.bind(i, tooltipMsg);
  return i;
}

export function extract(x: any): any {
  if (x == null)
    return null;
  if (x instanceof Widget)
    return x.root;
  if (typeof x.root !== 'undefined')
    return x.root;
  if (typeof x.wrap === 'function') {
    if (x.length === 1)
      return x[0];
    else
      return span(x.get());
  }  // jquery (cash) object
  return x;
}

/** Renders object to html element.
 * @param {object} x
 * @returns {HTMLElement} */
export function render(x: any): HTMLElement {
  x = extract(x);
  // @ts-ignore
  if (api.grok_UI_Render == null)
    return x;
  return api.grok_UI_Render(x);
}

/** Renders a table cell to html element, taking into account grid's formatting and color-coding.
 * @param {Cell} tableCell
 * @returns {HTMLElement} */
// export function renderCell(tableCell) {
//
// }

/** Renders a table cell to html element,
 * Takes into account grid's formatting and color-coding.
 * @param {Row} table
 * @param {string[]} columnNames
 * @returns {HTMLElement} */
// export function renderForm(table, columnNames = null) {
//   return null;
// }

/** Renders object to html card.
 * @param {object} x
 * @returns {HTMLElement}. */
export function renderCard(x: object): HTMLElement {
  return api.grok_UI_RenderCard(x);
}

/** Renders span
 * @param {object[]} x
 * @returns {HTMLElement}. */
export function span(x: object[]): HTMLElement {
  return api.grok_UI_Span(x);
}

export function renderInline(x: HTMLElement): HTMLElement {
  let handler = ObjectHandler.forEntity(x);
  return handler == null ? render(x) : handler.renderMarkup(x);
}


/** Renders inline text, calling [renderMarkup] for each non-HTMLElement
 * @param {object[]} objects
 * @returns {HTMLElement}. */
export function inlineText(objects: any[]): HTMLElement {
  return span(objects.map((item) => renderInline(item)));
}

/**
 * @param {object[]} children
 * @param {string | ElementOptions} options
 * @returns {HTMLDivElement}
 * */
export function div(children: HTMLElement[] = [], options: string | ElementOptions | null = null): HTMLDivElement {
  if (!Array.isArray(children))
    children = [children];
  let d = document.createElement('div');
  if (children != null) {
    $(d).append(children.map(render));
  }
  _options(d, options);
  $(d).addClass('ui-div');
  return d;
}

/** Div flex-box container that positions child elements vertically.
 * @param {object[]} items
 * @param {string | ElementOptions} options
 * @returns {HTMLDivElement} */
export function divV(items: HTMLElement[], options: string | ElementOptions | null = null): HTMLDivElement {
  return <HTMLDivElement>_options(api.grok_UI_DivV(items == null ? null : items.map(render), null), options);
}

/** Div flex-box container that positions child elements horizontally.
 * @param {object[]} items
 * @param {string | ElementOptions} options
 * @returns {HTMLDivElement} */
export function divH(items: HTMLElement[], options: string | ElementOptions | null = null): HTMLDivElement {
  return <HTMLDivElement>_options(api.grok_UI_DivH(items == null ? null : items.map(render), null), options);
}

/** Renders content as a card. */
export function card(content: HTMLElement): HTMLDivElement {
  return div([content], 'd4-item-card');
}

export function loader(): any {
  return api.grok_UI_Loader();
}

export function setUpdateIndicator(element: HTMLElement, updating: boolean = true): void {
  return api.grok_UI_SetUpdateIndicator(element, updating);
}

/**
 * Creates a button with the specified text, click handler, and tooltip
 * @param {string | Element | Array<string | Element>} content
 * @param {Function} handler
 * @param {string} tooltip
 * @returns {HTMLButtonElement}
 * */
export function button(content: string | Element | (string | Element)[], handler: Function, tooltip: string | null = null): HTMLButtonElement {
  return api.grok_UI_Button(content, handler, tooltip);
}

export function bigButton(text: string, handler: Function, tooltip: string | null = null): HTMLButtonElement {
  return api.grok_UI_BigButton(text, handler, tooltip);
}

/**
 * Creates a combo popup with the specified icons and items
 * @param {string | HTMLElement} caption
 * @param {Array<string>} items
 * @param {Function} handler (item) => {...}
 * @param {Function} renderer (item) => {...}
 * @returns {HTMLElement}
 * */
export function comboPopup(caption: string | HTMLElement, items: string[], handler: (item: any) => void, renderer: ((item: any) => HTMLElement) | null = null): HTMLElement {
  return api.grok_UI_ComboPopup(caption, items, handler, renderer !== null ? (item: any) => renderer(toJs(item)) : null);
}

/**
 * Creates a combo popup with the specified icons and items
 * @param {string | HTMLElement} caption
 * @param {Map<string>} items
 * @returns {HTMLElement}
 * */
export function comboPopupItems(caption: string | HTMLElement, items: { [key: string]: Function }): HTMLElement {
  return api.grok_UI_ComboPopup(caption, Object.keys(items), (key: string) => items[key]());
}

/** Creates a visual table based on [map]. */
export function tableFromMap(map: { [key: string]: any }): HTMLTableElement {
  return api.grok_UI_TableFromMap(map);
}

/** Creates a visual table based on [items] and [renderer]. */
export function table(items: [], renderer: ((item: any, ind: number) => any) | null, columnNames: string[] | null = null): HTMLTableElement {
  return api.grok_HtmlTable(items, renderer !== null ? (object: any, ind: number) => renderer(toJs(object), ind) : null, columnNames).root;
}

/** Waits for Future<Element> function to complete and collect its result.*/
export function wait(getElement: () => HTMLElement): any {
  return toJs(api.grok_UI_Wait(getElement));
}

/** Creates a visual element representing list of [items]. */
export function list(items: any[]): HTMLElement[] {
  return api.grok_UI_List(Array.from(items).map(toDart));
}

/** Creates a [Dialog].
 * @returns {Dialog} */
export function dialog(title: string = ''): Dialog {
  return Dialog.create(title);
}

/** Binds [item] with the [element]. It enables selecting it as a current object, drag-and-drop,
 * tooltip, and popup menu.
 * @param item
 * @param {Element} element
 * @returns {Element}. */
export function bind(item: any, element: Element): Element {
  return api.grok_UI_Bind(item, element);
}

/** Shows popup with the [element] near the [anchor].
 * tooltip, and popup menu.
 * @param {Element} element
 * @param {Element} anchor
 * @param {boolean} vertical
 * @returns {Element}. */
export function showPopup(element: Element, anchor: Element, vertical: boolean = false): Element {
  return api.grok_UI_ShowPopup(element, anchor, vertical);
}

export function rangeSlider(minRange: number, maxRange: number, min: number, max: number): RangeSlider {
  let rs = RangeSlider.create();
  rs.setValues(minRange, maxRange, min, max);
  return rs;
}

/**
 * Creates a virtual list widget
 * @param {number} length - number of elements
 * @param {Function} renderer
 * @param {boolean} verticalScroll - vertical or horizontal scrolling
 * @param {number} maxColumns - maximum number of items on the non-scrolling axis
 * @returns {VirtualView}
 */
export function virtualView(length: number, renderer: Function, verticalScroll: boolean = true, maxColumns: number = 1000): VirtualView {
  let view = VirtualView.create(verticalScroll, maxColumns);
  view.setData(length, renderer);
  return view;
}

export function popupMenu(items: any): void {
  function populate(menu: Menu, item: { [key: string]: any; }) {
    for (let key of Object.keys(item)) {
      let value = item[key];
      if (value instanceof Function)
        menu.item(key, value);
      else
        populate(menu.group(key), value);
    }
  }

  let menu = Menu.popup();
  populate(menu, items);
  menu.show();
}

export function inputs(inputs: InputBase[], options: {}) {
  return form(inputs, options);
}

/** Creates new nodes tree
 * @returns {TreeViewNode} */
export function tree(): TreeViewNode {
  return TreeViewNode.tree();
}

export function intInput(name: string, value: number, onValueChanged: Function | null = null): InputBase {
  return new InputBase(api.grok_IntInput(name, value), onValueChanged);
}

export function choiceInput<T>(name: string, selected: T, items: T[], onValueChanged: Function | null = null): InputBase {
  return new InputBase(api.grok_ChoiceInput(name, selected, items), onValueChanged);
}

export function multiChoiceInput<T>(name: string, value: T, items: T[], onValueChanged: Function | null = null): InputBase {
  return new InputBase(api.grok_MultiChoiceInput(name, value, items), onValueChanged);
}

export function stringInput(name: string, value: string, onValueChanged: Function | null = null): InputBase {
  return new InputBase(api.grok_StringInput(name, value), onValueChanged);
}

export function floatInput(name: string, value: number, onValueChanged: Function | null = null): InputBase {
  return new InputBase(api.grok_FloatInput(name, value), onValueChanged);
}

export function dateInput(name: string, value: DateTime, onValueChanged: Function | null = null): InputBase {
  return new InputBase(api.grok_DateInput(name, value.d), onValueChanged);
}

export function boolInput(name: string, value: boolean, onValueChanged: Function | null = null): InputBase {
  return new InputBase(api.grok_BoolInput(name, value), onValueChanged);
}

export function moleculeInput(name: string, value: string, onValueChanged: Function | null = null): InputBase {
  return new InputBase(api.grok_MoleculeInput(name, value), onValueChanged);
}

export function columnInput(name: string, table: DataFrame, value: Column, onValueChanged: Function | null = null): InputBase {
  return new InputBase(api.grok_ColumnInput(name, table.d, value.d), onValueChanged);
}

export function columnsInput(name: string, table: DataFrame, onValueChanged: Function | null = null): InputBase {
  return new InputBase(api.grok_ColumnsInput(name, table.d), onValueChanged);
}

export function textInput(name: string, value: string, onValueChanged: Function | null = null): InputBase {
  return new InputBase(api.grok_TextInput(name, value), onValueChanged);
}

/**
 * @param {HTMLElement} element
 * @returns {rxjs.Observable} */
export function onSizeChanged(element: HTMLElement): rxjs.Observable<any> {

  if (_isDartium()) {
    // Polyfill for Dartium which does not have a native ResizeObserver
    return new rxjs.Observable(function(observer) {
      let width = element.clientWidth;
      let height = element.clientHeight;
      let interval = setInterval(() => {
        let newWidth = element.clientWidth;
        let newHeight = element.clientHeight;
        if (newWidth !== width || newHeight !== height) {
          width = newWidth;
          height = newHeight;
          observer.next(element);
        }
      }, 100);
      return () => clearInterval(interval);
    });
  }

  return rxjs.Observable.create(function (observer: { next: (arg0: ResizeObserverEntry) => void; }) {
    const resizeObserver = new ResizeObserver(observerEntries => {
      // trigger a new item on the stream when resizes happen
      for (const entry of observerEntries) {
        observer.next(entry);
      }
    });

    // start listening for resize events
    resizeObserver.observe(element);

    // cancel resize observer on cancellation
    return () => resizeObserver.disconnect();
  });

}

/** UI Tools **/
export class tools {

  static handleResize(element: HTMLElement, onChanged: (width: number, height: number) => void): () => void {
    let width = element.clientWidth;
    let height = element.clientHeight;
    let interval = setInterval(() => {
      let newWidth = element.clientWidth;
      let newHeight = element.clientHeight;
      if (newWidth !== width || newHeight !== height) {
        width = newWidth;
        height = newHeight;
        onChanged(width, height);
      }
    }, 100);
    return () => clearInterval(interval);
  }
}

/** Represents a tooltip. */
export class Tooltip {

  /** Hides the tooltip. */
  hide(): void {
    api.grok_Tooltip_Hide();
  }

  /** Associated the specified visual element with the corresponding item. */
  bind(element: HTMLElement, tooltip: string): HTMLElement {
    api.grok_Tooltip_SetOn(element, tooltip);
    return element;
  }

  /** Shows the tooltip at the specified position
   * @param {HTMLElement | string} content
   * @param {number} x
   * @param {number} y */
  show(content: HTMLElement | string, x: number, y: number): void {
    api.grok_Tooltip_Show(content, x, y);
  }

  showRowGroup(dataFrame: DataFrame, indexPredicate: IndexPredicate, x: number, y: number): void {
    api.grok_Tooltip_ShowRowGroup(dataFrame.d, indexPredicate, x, y);
  }

  /** Returns a tooltip element.
   * @returns {HTMLElement} */
  get root(): HTMLElement {
    return api.grok_Tooltip_Get_Root();
  }

  /** @returns {boolean} isVisible */
  get isVisible(): boolean {
    return api.grok_Tooltip_Is_Visible();
  }

  /** @returns {Observable} */
  get onTooltipRequest(): rxjs.Observable<any> {
    return __obs('d4-tooltip-request');
  }

  /** @returns {Observable} */
  get onTooltipShown(): rxjs.Observable<any> {
    return __obs('d4-tooltip-shown');
  }

  /** @returns {Observable} */
  get onTooltipClosed(): rxjs.Observable<any> {
    return __obs('d4-tooltip-closed');
  }
}

export let tooltip = new Tooltip();

export class ObjectHandlerResolutionArgs {
  object: object;
  context: object | null;
  handler: null;

  constructor(object: object, context: object | null) {
    this.object = object;
    this.context = context;
    this.handler = null;
  }
}

let _objectHandlerSubject = new rxjs.Subject();

/**
 * Override the corresponding methods, and {@link register} an instance to
 * let Datagrok know how to handle objects of the specified type.
 *
 * {@link isApplicable} is used to associate an object with the handler.
 * When handling an object x, the platform uses the first registered handler that
 * claims that it is applicable to that object.
 *
 * TODO: search, destructuring to properties
 *
 * Samples: {@link https://public.datagrok.ai/js/samples/ui/meta/meta}
 * */
export class ObjectHandler {

  /** Type of the object that this meta handles. */
  get type(): string {
    throw 'Not defined.';
  }

  /**
   * Override this method to check whether this meta class should handle the specified object.
   * @param x - specified object.
   * @returns {boolean}
   * */
  isApplicable(x: any): boolean {
    throw 'Not defined.';
  }

  /** String representation of the [item], by default item.toString().
   * @param x - item
   * @returns {string} */
  getCaption(x: any): string {
    return `${x}`;
  }

  /** @returns {CanvasRenderer} */
  getCanvasRenderer(): CanvasRenderer | null {
    return null;
  }

  /** @returns {GridCellRenderer} */
  getGridCellRenderer(): GridCellRenderer | null {
    return null;
  }

  /** Renders icon for the item.
   * @param x - item
   * @returns {Element} */
  renderIcon(x: any, context: any = null): HTMLDivElement {
    return divText(this.getCaption(x));
  }

  /** Renders markup for the item.
   * @param x - item
   * @returns {Element} */
  renderMarkup(x: any, context: any = null): HTMLDivElement {
    return divText(this.getCaption(x));
  }

  /** Renders tooltip for the item.
   * @param x - item
   * @returns {Element} */
  renderTooltip(x: any, context: any = null): HTMLDivElement {
    return divText(this.getCaption(x));
  }

  /** Renders card div for the item.
   * @param x - item
   * @returns {Element} */
  renderCard(x: any, context: any = null): HTMLDivElement {
    return divText(this.getCaption(x));
  }

  /** Renders properties list for the item.
   * @param x - item
   * @returns {Element} */
  renderProperties(x: any, context: any = null): HTMLDivElement {
    return divText(this.getCaption(x));
  }

  /** Renders view for the item.
   * @param x - item
   * @returns {Element} */
  renderView(x: any, context: any = null): HTMLDivElement {
    return this.renderProperties(x);
  }

  /** Gets called once upon the registration of meta export class. */
  init(): void { }

  /** Registers entity handler.
   * @param {ObjectHandler} meta */
  static register(meta: ObjectHandler): void {
    api.grok_Meta_Register(meta);
    let cellRenderer = meta.getGridCellRenderer();
    if (cellRenderer != null)
      GridCellRenderer.register(cellRenderer);
  }

  /** @returns {ObjectHandler[]} */
  static list(): ObjectHandler[] {
    return toJs(api.grok_Meta_List());
  }

  static onResolve(observer: rxjs.PartialObserver<unknown> | undefined): rxjs.Subscription {
    return _objectHandlerSubject.subscribe(observer);
  }

  /**
   * @param {Object} object
   * @param {Object} context
   * @returns {ObjectHandler}
   * */
  static forEntity(object: object, context: object | null = null): ObjectHandler | null {
    let args = new ObjectHandlerResolutionArgs(object, context);
    _objectHandlerSubject.next(args);
    if (args.handler !== null)
      return args.handler;

    return toJs(api.grok_Meta_ForEntity(toDart(object)));
  }

  /**
   * Registers a function that takes applicable objects as the only argument.
   * It will be suggested to run in the context menu for that object, and
   * also in the "Actions" pane on the property panel.
   *
   * Samples: {@link https://public.datagrok.ai/js/samples/ui/docking/docking}
   *
   * @param {string} name - function name
   * @param run - a function that takes exactly one parameter
   * */
  registerParamFunc(name: string, run: (param: any) => any): void {
    // @ts-ignore
    new Functions().registerParamFunc(name, this.type, run, this.isApplicable);
  }
}

export class EntityMetaDartProxy extends ObjectHandler {
  d: any;

  constructor(d: any) {
    super();
    this.d = d;
  }

  get type(): string { return api.grok_Meta_Get_Type(this.d); }
  isApplicable(x: any): boolean { return api.grok_Meta_IsApplicable(this.d, toDart(x)); }
  getCaption(x: any): string { return api.grok_Meta_Get_Name(this.d, x); }

  renderIcon(x: any, context = null): HTMLDivElement { return api.grok_Meta_RenderIcon(x); }
  renderMarkup(x: any, context = null): HTMLDivElement { return api.grok_Meta_RenderMarkup(x); }
  renderTooltip(x: any, context = null): HTMLDivElement { return api.grok_Meta_RenderTooltip(x); }
  renderCard(x: any, context = null): HTMLDivElement { return api.grok_Meta_RenderCard(x); }
  renderProperties(x: any, context = null): HTMLDivElement { return api.grok_Meta_RenderProperties(x); }
  renderView(x: any, context = null): HTMLDivElement { return api.grok_Meta_RenderProperties(x); }
}

export function box(item: Widget | InputBase | HTMLElement | null, options: ElementOptions | null = null): HTMLDivElement {
  if (item instanceof Widget) {
    item = item.root;
  }
  if (item instanceof InputBase) {
    console.log('inputbase');
    item = item.root;
  }
  if (item != null && $(item).hasClass('ui-box')) {
    return item as HTMLDivElement;
  }
  let c = document.createElement('div');
  _options(c, options);
  $(c).append(item).addClass('ui-box');
  return c;
}

export function boxFixed(item: Widget | InputBase | HTMLElement | null, options: ElementOptions | null = null): HTMLDivElement {
  let c = box(item, options);
  $(c).addClass('ui-box-fixed');
  return c;
}

/** Div flex-box container that positions child elements vertically.
 *  @param {object[]} items */
export function splitV(items: HTMLElement[], options: ElementOptions | null = null): HTMLDivElement {
  let b = box(null, options);
  $(b).addClass('ui-split-v').append(items.map(item => box(item)));
  return b;
}

/** Div flex-box container that positions child elements horizontally.
 *  @param {object[]} items */
export function splitH(items: HTMLElement[], options: ElementOptions | null = null): HTMLDivElement {
  let b = box(null, options);
  $(b).addClass('ui-split-h').append(items.map(item => box(item)));
  return b;
}

export function block(items: HTMLElement[], options: ElementOptions | null = null): HTMLDivElement {
  let c = div(items, options);
  $(c).addClass('ui-block');
  return c;
}

export function block75(items: HTMLElement[], options: ElementOptions | null = null) {
  return $(block(items, options)).addClass('ui-block-75')[0];
}
export function block25(items: HTMLElement[], options: ElementOptions | null = null) {
  return $(block(items, options)).addClass('ui-block-25')[0];
}
export function block50(items: HTMLElement[], options: ElementOptions | null = null) {
  return $(block(items, options)).addClass('ui-block-50')[0];
}

export function p(text: string, options = null): HTMLParagraphElement {
  let c = document.createElement('p');
  c.textContent = text;
  $(c).addClass('ui-p');
  _options(c, options);
  return c;
}

export function panel(items: HTMLElement[], options?: string | ElementOptions) {
  return $(div(items, options)).addClass('ui-panel')[0];
}

export function label(text: string | null, options: {} | null = null): HTMLLabelElement {
  let c = document.createElement('label');
  c.textContent = text;
  $(c).addClass('ui-label');
  _options(c, options);
  return c;
}

export function form(children: InputBase[] = [], options: {} | null = null): HTMLDivElement {
  let d = document.createElement('div');
  if (children != null) {
    $(d).append(children.map(render));
  }
  _options(d, options);
  $(d).addClass('ui-form');
  return d;
}

export function narrowForm(children: InputBase[] = [], options: {} | null = null): HTMLDivElement {
  let d = form(children, options);
  $(d).addClass('ui-form-condensed');
  return d;
}

export function buttonsInput(children: HTMLButtonElement[] = []): HTMLDivElement {
  if (!Array.isArray(children))
    children = [children];

  let d = document.createElement('div');
  let l = document.createElement('label');
  let e = document.createElement('div');
  $(e).addClass('ui-input-editor');
  if (children != null)
    $(e).append(children.map(render));
  l.textContent = ' ';
  $(l).addClass('ui-label ui-input-label');
  $(d).addClass('ui-input-root ui-input-buttons').append(l).append(e);
  return d;
}

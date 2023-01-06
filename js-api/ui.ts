/**
 * Routines for building UI
 * @module ui
 **/

import {ElementOptions, IndexPredicate} from './src/const';
import {Viewer} from './src/viewer';
import {VirtualView} from './src/views/view';
import {
  Accordion,
  Dialog,
  InputBase,
  Menu,
  TabControl,
  TreeViewGroup,
  Widget,
  RangeSlider,
  RangeSliderStyle,
  FilesWidget,
  DateInput,
  fileShares,
  SliderOptions
} from './src/widgets';
import {toDart, toJs} from './src/wrappers';
import {Functions} from './src/functions';
import $ from 'cash-dom';
import {__obs, StreamSubscription} from './src/events';
import {_isDartium, _options} from './src/utils';
import * as rxjs from 'rxjs';
import { CanvasRenderer, GridCellRenderer, SemanticValue } from './src/grid';
import {Entity, Property} from './src/entities';
import { Column, DataFrame } from './src/dataframe';
import dayjs from "dayjs";


let api = <any>window;
declare let grok: any;

/**
 * Creates an instance of the element for the specified tag, and optionally assigns it a CSS class.
 * @param {string} tagName The name of an element.
 * @param {string | null} className
 * @returns {HTMLElement}
 */
export function element(tagName: string, className: string | null = null): HTMLElement & any {
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

/** Removes all child nodes, and disposes widgets that these nodes
 *  might contain. See also [remove]. */
export function empty(e: HTMLElement): HTMLElement {
  while (e.firstChild) {
    api.grok_Widget_Kill(e.firstChild);
    e.removeChild(e.firstChild);
  }
  return e;
}

export function setClass(e: HTMLElement, classes: string, flag: boolean) {
  if (flag)
    $(e).addClass(classes);
  else
    $(e).removeClass(classes);
}

/** Removes the [element] from the DOM, and disposes any widgets that
 *  the [element] might contain. See also [empty]. */
export function remove(element: HTMLElement): void {
  api.grok_Widget_Kill(element);
  element.remove();
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
export function h1(s: string | Element, options: string | ElementOptions | null = null): HTMLHeadingElement {
  let x = element('h1');
  if (typeof s === 'string')
    x.innerText = s;
  else
    $(x).append(s);
  return _options(x, options) as HTMLHeadingElement;
}

/** @returns {HTMLHeadingElement} */
export function h2(s: string | Element, options: string | ElementOptions | null = null): HTMLHeadingElement {
  let x = element('h2');
  if (typeof s === 'string')
    x.innerText = s;
  else
    $(x).append(s);
  return _options(x, options) as HTMLHeadingElement;
}

/** @returns {HTMLHeadingElement} */
export function h3(s: string | Element, options: string | ElementOptions | null = null): HTMLHeadingElement {
  let x = element('h3');
  if (typeof s === 'string')
    x.innerText = s;
  else
    $(x).append(s);
  return _options(x, options) as HTMLHeadingElement;
}

/** Creates an accordion with dynamically populated panes.
 * Sample: {@link https://public.datagrok.ai/js/samples/ui/ui}
 * @returns {Accordion} */
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
export function divText(text: string, options: string | ElementOptions | any | null = null): HTMLDivElement {
  let e = element('div');
  e.innerText = text;
  _options(e, options);
  return e as HTMLDivElement;
}

export function tags(entity: Entity): HTMLElement {
  return api.grok_UI_Tags(entity.dart);
}

export function markdown(text: string): HTMLElement {
  return api.grok_UI_Markdown(text);
}

/** Returns a font-awesome icon with the specified name, handler, and tooltip.
 * Sample: {@link https://public.datagrok.ai/js/samples/ui/icons}
 * @param {string} name - icon name (omit the "fa-" prefix)
 * @param {Function} handler
 * @param {String} tooltipMsg
 * @returns {HTMLElement} */
export function iconFA(name: string, handler: ((this: HTMLElement, ev: MouseEvent) => any) | null = null, tooltipMsg: string | null = null): HTMLElement {
  let i = element('i');
  i.classList.add('grok-icon');
  i.classList.add('fal');
  i.classList.add(`fa-${name}`);
  if (handler !== null)
    i.addEventListener('click', handler);
  if (tooltipMsg !== null)
    tooltip.bind(i, tooltipMsg);
  return i;
}

export function iconImage(name: string, path: string, handler: ((this: HTMLElement, ev: MouseEvent) => any) | null = null, tooltipMsg: string | null = null): HTMLElement {
  let i = element('i');
  i.classList.add('grok-icon');
  i.classList.add('image-icon');
  if (!path.startsWith('http') && !path.startsWith('/') && !path.startsWith('data:'))
    path = '/images/$path';
  i.style.backgroundImage = `url(${path})`;
  if (handler !== null)
    i.addEventListener('click', handler);
  if (tooltipMsg !== null)
    tooltip.bind(i, tooltipMsg);
  return i;
}

export function iconSvg(name: string, handler: ((this: HTMLElement, ev: MouseEvent) => any) | null = null, tooltipMsg: string | null = null): HTMLElement {
  let i = element('i');
  i.classList.add('svg-icon');
  i.classList.add('grok-icon');
  i.classList.add(`svg-${name}`);
  if (handler !== null)
    i.addEventListener('click', handler);
  if (tooltipMsg !== null)
    tooltip.bind(i, tooltipMsg);
  return i;
}

export function extractRoot(x: any): HTMLElement | null {
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

//function applyToNode(HTMLElement element)

/** Renders object to html element. */
export function render(x: any, options?: ElementOptions): HTMLElement {
  function renderImpl() {
    x = extractRoot(x);
    if (x == null)
      return div();
    if (api.grok_UI_Render == null)
      return x;
    return api.grok_UI_Render(x);
  }
  return _options(renderImpl(), options);
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

/** Renders multiple objects as a span */
export function span(x: any[], options: string | ElementOptions | null = null): HTMLElement {
  return _options(api.grok_UI_Span(x), options);
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
export function div(children: any[] | string | HTMLElement = [], options: string | ElementOptions | null = null): HTMLDivElement {
  if (!Array.isArray(children))
    children = [children];
  let d = document.createElement('div');
  if (children != null) {
    $(d).append(children.map(x => render(x)));
  }
  _options(d, options);
  $(d).addClass('ui-div');
  return d;
}

export function info(children: HTMLElement[] | HTMLElement | string, header: string | null = null, reopenable: boolean = true): HTMLDivElement {
  let root: HTMLDivElement | null;
  let divContent: HTMLElement[] = [];
  let divActual: HTMLDivElement | null = null;
  let show = iconFA('lightbulb-exclamation', () => {
    if (divActual) divActual.style.display = 'block';
    close.style.display = 'block';
    show.style.display = 'none';
  });
  let close = iconFA('times', () => {
    if (divActual) divActual.style.display = 'none';
    close.style.display = 'none';
    show.style.display = reopenable ? 'block' : 'none';
  });
  if (header !== null && header !== undefined) {
    divContent.push(h1(header));
  }
  if (!Array.isArray(children)) {
    if (children === null || typeof children === 'string') {
      children = [divText(<string>children)];
    } else {
      children = [children];
    }
  }
  divContent = divContent.concat(children);
  divActual = div(divContent, 'grok-info-bar');
  root = div([show, close, divActual], 'grok-info-bar-container');
  return root;
}

/** Div flex-box container that positions child elements vertically.
 * @param {object[]} items
 * @param {string | ElementOptions} options
 * @returns {HTMLDivElement} */
export function divV(items: any[], options: string | ElementOptions | null = null): HTMLDivElement {
  return <HTMLDivElement>_options(api.grok_UI_DivV(items == null ? null : items.map(x => render(x)), 'ui-div'), options);
}

/** Div flex-box container that positions child elements horizontally.
 * @param {object[]} items
 * @param {string | ElementOptions} options
 * @returns {HTMLDivElement} */
export function divH(items: HTMLElement[], options: string | ElementOptions | null = null): HTMLDivElement {
  return <HTMLDivElement>_options(api.grok_UI_DivH(items == null ? null : items.map(x => render(x)), 'ui-div'), options);
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
export function comboPopup(caption: string | HTMLElement, items: string[], handler: (item: any) => void, renderer?: ((item: any) => HTMLElement) | null): HTMLElement {
  return api.grok_UI_ComboPopup(caption, items, handler, renderer ? (item: any) => renderer(toJs(item)) : null);
}

/**
 * Creates a combo popup with the specified icons and items
 * @param {string | HTMLElement} caption
 * @param {Map<string>} items
 * @returns {HTMLElement}
 * */
export function comboPopupItems(caption: string | HTMLElement, items: { [key: string]: Function }): HTMLElement {
  return api.grok_UI_ComboPopup(caption, Object.keys(items), (key: string) => items[key](), null);
}

/** Creates an html table based on [map]. */
export function tableFromMap(map: { [key: string]: any }): HTMLTableElement {
  return api.grok_UI_TableFromMap(map);
}

/** Creates an editable html table for the specified items (rows) and properties (columns). */
export function tableFromProperties(items: any[], properties: Property[]) {
  return table(items,
    (item, _) => properties.map((p) => InputBase.forProperty(p, item).input),
    properties.map((p) => p.name));
}

/** Creates a visual table based on [items] and [renderer]. */
export function table<T>(items: T[], renderer: ((item: T, ind: number) => any) | null, columnNames: string[] | null = null): HTMLTableElement {
  return toJs(api.grok_HtmlTable(items, renderer !== null ? (object: any, ind: number) => renderer(toJs(object), ind) : null, columnNames)).root;
}

/** Waits for Future<Element> function to complete and collect its result.*/
export function wait(getElement: () => Promise<HTMLElement>): any {
  return toJs(api.grok_UI_Wait(getElement));
}

/** Waits for Future<Element> function to complete and collect its result as a ui.box.*/
export function waitBox(getElement: () => Promise<HTMLElement>): any {
  return toJs(api.grok_UI_WaitBox(getElement));
}

/** Creates a visual element representing list of [items]. */
export function list(items: any[], options?: {processNode?: (node: HTMLElement) => void}): HTMLElement {
  const host: HTMLElement = api.grok_UI_List(Array.from(items).map(toDart));
  if (options?.processNode != null)
    for (const c of Array.from(host.children))
      options.processNode(c as HTMLElement);
  return host;
}

export function iframe(options?: {src?: string, width?: string, height?: string}) {
  let frame = element('iframe') as HTMLIFrameElement;
  if (options?.src != null)
    frame.src = options?.src;
  if (options?.width != null) {
    frame.width = options.width;
    frame.style.width = options.width;
  }
  if (options?.height != null) {
    frame.height = options.height;
    frame.style.height = options.height;
  }

  return frame;
}

function _link(element: HTMLElement, target: string | Function, tooltipMsg?: string): void {
  if (!tooltipMsg && typeof target === 'string')
    tooltipMsg = 'Open in new tab';

  element.style.cursor = 'pointer';
  element.addEventListener('click', (_) => {
    if (target instanceof Function)
      target();
    else
      window.open(target, '_blank');
  });

  tooltip.bind(element, tooltipMsg);
}

export function image(src: string, width: number, height: number, options?: {target?: string | Function, tooltipMsg?: string}) {
  let image = element('div') as HTMLDivElement;
  image.classList.add('ui-image');

  image.style.backgroundImage = `url(${src})`;
  image.style.width = `${width}px`;
  image.style.height = `${height}px`;

  if (options?.target)
    _link(image, options?.target, options?.tooltipMsg);

  return image;
}

/** Creates an <a> element. */
export function link(
    text: string,
    target: string | Function | object,
    tooltipMsg?: string,
    options: string | ElementOptions | null = { classes: 'd4-link-external' }): HTMLAnchorElement {
  let link = element('a') as HTMLAnchorElement;
  link.classList.add('ui-link');
  link.innerText = text;
  if (target instanceof Function) {
    link.addEventListener('click', (_) => target());
  }
  else if (typeof target === 'string') {
    link.href = target;
    link.target = '_blank';
  }
  else {
    bind(target, link);
  }
  tooltip.bind(link, tooltipMsg);
  _options(link, options);
  return link;
}

/** Creates a [Dialog]. */
export function dialog(options?: { title?: string, helpUrl?: string } | string): Dialog {
  return Dialog.create(options);
}

/** Binds [item] with the [element]. It enables selecting it as a current object, drag-and-drop,
 * tooltip, and popup menu.*/
export function bind(item: any, element: HTMLElement, options?: {contextMenu: boolean}): HTMLElement {
  return api.grok_UI_Bind(toDart(item), element, options?.contextMenu);
}

/** Shows popup with the [element] near the [anchor].
 * tooltip, and popup menu.
 * @param {Element} element
 * @param {Element} anchor
 * @param {boolean} vertical
 * @returns {Element}. */
export function showPopup(element: HTMLElement, anchor: HTMLElement, vertical: boolean = false): Element {
  return api.grok_UI_ShowPopup(element, anchor, vertical);
}

export function rangeSlider(minRange: number, maxRange: number, min: number, max: number, vertical: boolean = false, style: RangeSliderStyle | SliderOptions = 'barbell'): RangeSlider {
  let rs = RangeSlider.create(vertical, style);
  rs.setValues(minRange, maxRange, min, max);
  return rs;
}

/**
 * Creates a virtual list widget.
 * Sample: {@link https://public.datagrok.ai/js/samples/ui/virtual-view}
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

export function makeDraggable<T>(e: Element,
    options?: {
      allowCopy?: () => boolean,
      check?: () => boolean,
      getDragObject?: () => any,
      getDragCaption?: () => any,
      dragObjectType?: string,
      getDragHint?: () => string,
      getDragContext?: () => any,
      onDragStart?: (me: MouseEvent) => boolean,
      onDragEnd?: () => void
    }): Element {
  return api.grok_UI_MakeDraggable(e,
      options?.allowCopy,
      options?.check,
      options?.getDragObject ? () => toDart(options.getDragObject!()) : null,
      options?.getDragCaption ? () => toDart(options.getDragCaption!()) : null,
      options?.dragObjectType,
      options?.getDragHint,
      options?.getDragContext ? () => toDart(options.getDragContext!()) : null,
      options?.onDragStart,
      options?.onDragEnd
  );
}

export function makeDroppable<T>(e: Element,
    options?: {
      acceptDrop?: (dragObject: T) => boolean,
      doDrop?: (dragObject: T, copying: boolean) => void
    }): void {
  api.grok_UI_MakeDroppable(e,
      options?.acceptDrop ? (dragObject: T) => options.acceptDrop!(toJs(dragObject)) : null,
      options?.doDrop ? (dragObject: T, copying: boolean) => options.doDrop!(toJs(dragObject), toJs(copying)) : null,
  );
}

export function bindInputs(inputs: InputBase[]): StreamSubscription[] {
  let s: StreamSubscription[] = [];
  inputs.map((i) => {
    inputs.map((j) => {
      if (j != i)
        s.push(j.onChanged(() => {
          i.notify = false;
          i.stringValue = j.stringValue;
          i.notify = true;
        }));
    })
  });
  return s;
}

export function inputs(inputs: Iterable<InputBase>, options: any = null) {
  return form([...inputs], options);
}

/** Creates new nodes tree.
 * Sample: {@link https://public.datagrok.ai/js/samples/ui/tree-view}
 * @returns {TreeViewGroup} */
export function tree(): TreeViewGroup {
  return TreeViewGroup.tree();
}


// Will come back later (-Andrew)
//
// export interface InputOptions {
//   value?: any;
//   onValueChanged?: (v: any) => void;
// }
//
//
export namespace input {

  /** Creates input for the specified property, and optionally binds it to the specified object */
  export function forProperty(property: Property, source: any = null): InputBase {
    return InputBase.forProperty(property, source);
  }

  /** Returns a form for the specified properties, bound to the specified object */
  export function form(source: any, props: Property[]): HTMLElement {
    return inputs(props.map((p) => forProperty(p, source)));
  }

  // export function bySemType(semType: string) {
  //
  // }

  // function _options(inputBase: InputBase, options?: InputOptions) {
  //   if (options == null)
  //     return;
  //   return inputBase;
  // }
  //
  // export function int(name: string, options?: InputOptions) {
  //   return _options(intInput(name, options?.value, options?.onValueChanged));
  // }
  //
  // export function choice(name: string, choices: string[], options?: InputOptions) {
  //   return _options(choiceInput(name, options?.value ?? choices[0], choices, options?.onValueChanged));
  // }
  //
  // export function multiChoice(name: string, options?: InputOptions) {
  //
  // }
}

export function inputsRow(name: string, inputs: InputBase[]): HTMLElement {
  let d = div([label(name, 'ui-label ui-input-label')], 'ui-input-root ui-input-row');
  d.appendChild(div(inputs));
  return d;
}

export function intInput(name: string, value: number | null, onValueChanged: Function | null = null): InputBase<number | null> {
  return new InputBase(api.grok_IntInput(name, value), onValueChanged);
}

export function sliderInput(name: string, value: number | null, min: number, max: number, onValueChanged: Function | null = null): InputBase<number | null> {
  return new InputBase(api.grok_SliderInput(name, value, min, max), onValueChanged);
}

export function choiceInput<T>(name: string, selected: T, items: T[], onValueChanged: Function | null = null): InputBase<T | null> {
  return new InputBase(api.grok_ChoiceInput(name, selected, items), onValueChanged);
}

export function multiChoiceInput<T>(name: string, value: T[], items: T[], onValueChanged: Function | null = null): InputBase<T[] | null> {
  return new InputBase(api.grok_MultiChoiceInput(name, value, items), onValueChanged);
}

export function stringInput(name: string, value: string, onValueChanged: Function | null = null, options: { icon?: string | HTMLElement, clearIcon?: boolean, escClears?: boolean, placeholder?: String } | null = null): InputBase<string> {
  return new InputBase(api.grok_StringInput(name, value, options), onValueChanged);
}

export function searchInput(name: string, value: string, onValueChanged: Function | null = null): InputBase<string> {
  return new InputBase(api.grok_SearchInput(name, value), onValueChanged);
}

export function floatInput(name: string, value: number | null, onValueChanged: Function | null = null): InputBase<number | null> {
  return new InputBase(api.grok_FloatInput(name, value), onValueChanged);
}

export function dateInput(name: string, value: dayjs.Dayjs, onValueChanged: Function | null = null): DateInput {
  return new DateInput(api.grok_DateInput(name, value.valueOf()), onValueChanged);
}

export function boolInput(name: string, value: boolean, onValueChanged: Function | null = null): InputBase<boolean | null> {
  return new InputBase(api.grok_BoolInput(name, value), onValueChanged);
}

export function switchInput(name: string, value: boolean, onValueChanged: Function | null = null): InputBase<boolean> {
  return new InputBase(api.grok_SwitchInput(name, value), onValueChanged);
}

export function moleculeInput(name: string, value: string, onValueChanged: Function | null = null): InputBase<string> {
  return new InputBase(api.grok_MoleculeInput(name, value), onValueChanged);
}

export function columnInput(name: string, table: DataFrame, value: Column | null, onValueChanged: Function | null = null): InputBase<Column | null> {
  return new InputBase(api.grok_ColumnInput(name, table.dart, value?.dart), onValueChanged);
}

export function columnsInput(name: string, table: DataFrame, onValueChanged: (columns: Column[]) => void,
                             options?: {available?: string[], checked?: string[]}): InputBase<Column[]> {
  return new InputBase(api.grok_ColumnsInput(name, table.dart, options?.available, options?.checked), onValueChanged);
}

export function tableInput(name: string, table: DataFrame | null, tables: DataFrame[] = grok.shell.tables, onValueChanged: Function | null = null): InputBase<DataFrame | null> {
  return new InputBase(api.grok_TableInput(name, table, tables), onValueChanged);
}

export function textInput(name: string, value: string, onValueChanged: Function | null = null): InputBase<string> {
  return new InputBase(api.grok_TextInput(name, value), onValueChanged);
}

export function colorInput(name: string, value: string, onValueChanged: Function | null = null): InputBase<string> {
  return new InputBase(api.grok_ColorInput(name, value), onValueChanged);
}

/**
 * Sample: {@link https://public.datagrok.ai/js/samples/ui/ui-events}
 * @param {HTMLElement} element
 * @returns {rxjs.Observable} */
export function onSizeChanged(element: HTMLElement): rxjs.Observable<any> {

  if (_isDartium()) {
    // Polyfill for Dartium which does not have a native ResizeObserver
    return new rxjs.Observable(function(observer: any) {
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

  static setHoverVisibility(host: HTMLElement, items: HTMLElement[]) {
    api.grok_Tools_SetHoverVisibility(host, items);
  }

  static createElementFromHtml(htmlString: string): HTMLDivElement {
    let host = div([]);
    host.innerHTML = htmlString.trim();
    return host;
  }

  static initFormulaAccelerators(textInput: InputBase, table: DataFrame): void {
    api.grok_UI_InitFormulaAccelerators(toDart(textInput), table.dart);
  }

  /** Waits until the specified element is in the DOM. */
  static waitForElementInDom(element: HTMLElement): Promise<HTMLElement> {
    if (_isDartium()) {
      return new Promise(resolve => {
        setInterval(function() {
          if (document.contains(element))
            return resolve(element);
        }, 100);
      });
    }

    return new Promise(resolve => {
      if (document.contains(element))
        return resolve(element);

      const observer = new MutationObserver(_ => {
        if (document.contains(element)) {
          resolve(element);
          observer.disconnect();
        }
      });

      observer.observe(document.body, {
        childList: true,
        subtree: true
      });
    });
  }
}

/** Represents a tooltip. */
export class Tooltip {

  /** Hides the tooltip. */
  hide(): void {
    api.grok_Tooltip_Hide();
  }

  /** Associated the specified visual element with the corresponding item. */
  bind(element: HTMLElement, tooltip?: string | null | (() => string | HTMLElement | null)): HTMLElement {
    if (tooltip != null)
      api.grok_Tooltip_SetOn(element, tooltip);
    return element;
  }

  /** Shows the tooltip at the specified position */
  show(content: HTMLElement | string, x: number, y: number): void {
    api.grok_Tooltip_Show(content, x, y);
  }

  showRowGroup(dataFrame: DataFrame, indexPredicate: IndexPredicate, x: number, y: number): void {
    api.grok_Tooltip_ShowRowGroup(dataFrame.dart, indexPredicate, x, y);
  }

  /** Returns a tooltip element. */
  get root(): HTMLElement {
    return api.grok_Tooltip_Get_Root();
  }

  get isVisible(): boolean { return api.grok_Tooltip_Is_Visible(); }

  get onTooltipRequest(): rxjs.Observable<any> { return __obs('d4-tooltip-request'); }
  get onTooltipShown(): rxjs.Observable<any> { return __obs('d4-tooltip-shown'); }
  get onTooltipClosed(): rxjs.Observable<any> { return __obs('d4-tooltip-closed'); }
}

export let tooltip = new Tooltip();

export class ObjectHandlerResolutionArgs {
  semValue: SemanticValue;
  context: object | null;
  handler: ObjectHandler | null;

  get value(): any { return this.semValue.value; }

  constructor(semValue: SemanticValue, context: object | null) {
    this.semValue = semValue;
    this.context = context;
    this.handler = null;
  }
}

let _objectHandlerSubject = new rxjs.Subject<ObjectHandlerResolutionArgs>();

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

  get name(): string { return `${this.type} handler`; }

  toString(): string { return this.name; }

  async getById(id: string): Promise<any> {
    return null;
  }

  async refresh(x: any): Promise<any> {
    return x;
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

  /** Renders icon for the item. */
  renderIcon(x: any, context: any = null): HTMLElement {
    return divText(this.getCaption(x));
  }

  /** Renders markup for the item. */
  renderMarkup(x: any, context: any = null): HTMLElement {
    return divText(this.getCaption(x));
  }

  /** Renders tooltip for the item. */
  renderTooltip(x: any, context: any = null): HTMLElement {
    return divText(this.getCaption(x));
  }

  /** Renders card div for the item. */
  renderCard(x: any, context: any = null): HTMLElement {
    return divText(this.getCaption(x));
  }

  /** Renders properties list for the item. */
  renderProperties(x: any, context: any = null): HTMLElement {
    return divText(this.getCaption(x));
  }

  /** Renders view for the item. */
  renderView(x: any, context: any = null): HTMLElement {
    return this.renderProperties(x);
  }

  /** Gets called once upon the registration of meta export class. */
  init(): void { }

  /** Registers entity handler. */
  static register(meta: ObjectHandler): void {
    api.grok_Meta_Register(meta);
    let cellRenderer = meta.getGridCellRenderer();
    if (cellRenderer != null)
      GridCellRenderer.register(cellRenderer);
  }

  static list(): ObjectHandler[] {
    return toJs(api.grok_Meta_List());
  }

  static onResolve(observer: rxjs.PartialObserver<ObjectHandlerResolutionArgs>): rxjs.Subscription {
    return _objectHandlerSubject.subscribe(observer);
  }

  static forEntity(object: object, context: object | null = null): ObjectHandler | null {
    let semValue = object instanceof SemanticValue ? object : SemanticValue.fromValueType(object, null);

    let args = new ObjectHandlerResolutionArgs(semValue, context);
    _objectHandlerSubject.next(args);
    if (args.handler !== null)
      return args.handler;

    return toJs(api.grok_Meta_ForEntity(toDart(object)));
  }

  static forSemType(semType: string): Promise<ObjectHandler[]> {
    return new Promise((resolve, reject) => api.grok_Meta_ForSemType(semType, (t: any) => resolve(toJs(t)), (e: any) => reject(e)));
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
  dart: any;

  constructor(d: any) {
    super();
    this.dart = d;
  }

  get type(): string { return api.grok_Meta_Get_Type(this.dart); }
  isApplicable(x: any): boolean { return api.grok_Meta_IsApplicable(this.dart, toDart(x)); }
  getCaption(x: any): string { return api.grok_Meta_Get_Name(this.dart, x); }

  renderIcon(x: any, context: any = null): HTMLDivElement { return api.grok_Meta_RenderIcon(x); }
  renderMarkup(x: any, context: any = null): HTMLDivElement { return api.grok_Meta_RenderMarkup(x); }
  renderTooltip(x: any, context: any = null): HTMLDivElement { return api.grok_Meta_RenderTooltip(x); }
  renderCard(x: any, context: any = null): HTMLDivElement { return api.grok_Meta_RenderCard(x); }
  renderProperties(x: any, context: any = null): HTMLDivElement { return api.grok_Meta_RenderProperties(x); }
  renderView(x: any, context: any = null): HTMLDivElement { return api.grok_Meta_RenderProperties(x); }
}

export function box(item: Widget | InputBase | HTMLElement | null = null, options: string | ElementOptions | null = null): HTMLDivElement {
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

/** Div flex-box container that positions child elements vertically. */
export function splitV(items: HTMLElement[], options: ElementOptions | null = null, resize: boolean | null = false): HTMLDivElement {
  let b = box(null, options);
  if (resize && items.length > 1){
    items.forEach((v, i) => {
      let divider = box();
      $(b).addClass('ui-split-v').append(box(v))
      if (i != items.length - 1) {
        $(b).append(divider)
        spliterResize(divider, items[i], items[i + 1])
      }
    })
  } else {
    $(b).addClass('ui-split-v').append(items.map(item => box(item)))
  }
  return b;
}

/** Div flex-box container that positions child elements horizontally. */
export function splitH(items: HTMLElement[], options: ElementOptions | null = null, resize: boolean | null = false): HTMLDivElement {
  let b = box(null, options);
  if (resize && items.length > 1) {
    items.forEach((v, i) => {
      let divider = box();
      $(b).addClass('ui-split-h').append(box(v))
      if (i != items.length - 1) {
        $(b).append(divider)
        spliterResize(divider, items[i], items[i + 1], true)
      }
    })
  } else {
    $(b).addClass('ui-split-h').append(items.map(item => box(item)))
  }
  return b;
}

function spliterResize(divider: HTMLElement, previousSibling: HTMLElement, nextSibling: HTMLElement, horizontal: boolean = false) {
  let md: any;
  divider.onmousedown = onMouseDown;

  if (horizontal) {
    divider.style.cssText = `
    background-image: url("data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' width='30' height='20'><path d='M 8 3 h 10 M 8 6 h 10 M 8 9 h 10' fill='none' stroke='%239497A0' stroke-width='1.25'/></svg>");
    max-width: 4px;
    cursor: col-resize;`
  }
  else {
    divider.style.cssText = `
    background-image: url("data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' width='30' height='20'><path d='M 2 5 v 10 M 5 5 v 10 M 8 5 v 10' fill='none' stroke='%239497A0' stroke-width='1.25'/></svg>");
    max-height: 4px;
    cursor: row-resize;`
  }

  divider.style.backgroundRepeat = 'no-repeat';
  divider.style.backgroundPosition = 'center';
  divider.style.backgroundColor = 'var(--grey-1)';

  function onMouseDown(e: any) {
    if (!nextSibling.classList.contains('ui-box'))
      nextSibling = nextSibling.parentElement!;
    if (!previousSibling.classList.contains('ui-box'))
      previousSibling = previousSibling.parentElement!;

    md = {
      e,
      offsetLeft: divider.offsetLeft,
      offsetTop: divider.offsetTop,
      topHeight: previousSibling.offsetHeight,
      bottomHeight: nextSibling.offsetHeight,
      leftWidth: previousSibling.offsetWidth,
      rightWidth: nextSibling.offsetWidth
    };
    divider.style.backgroundColor = 'var(--grey-2)';
    document.onmousemove = onMouseMove;
    document.onmouseup = () => {
      divider.style.backgroundColor = 'var(--grey-1)';
      document.onmousemove = document.onmouseup = null;
    }
  }

  function onMouseMove(e: any) {
    const delta = {x: e.clientX - md.e.clientX, y: e.clientY - md.e.clientY};
    if (horizontal) {
      delta.x = Math.min(Math.max(delta.x, -md.leftWidth), md.rightWidth);
      previousSibling.style.maxWidth = (md.leftWidth + delta.x) + "px";
      nextSibling.style.maxWidth = (md.rightWidth - delta.x) + "px";
    } else {
      delta.x = Math.min(Math.max(delta.y, -md.topHeight), md.bottomHeight);
      previousSibling.style.maxHeight = (md.topHeight + delta.y) + "px";
      nextSibling.style.maxHeight = (md.bottomHeight - delta.y) + "px";
    }
  }
}

export function block(items: HTMLElement[], options: string | ElementOptions | null = null): HTMLDivElement {
  let c = div(items, options);
  $(c).addClass('ui-block');
  return c;
}

export function block75(items: HTMLElement[], options: ElementOptions | null = null): HTMLDivElement {
  let c = block(items, options);
  $(c).addClass('ui-block-75');
  return c;
}

export function block25(items: HTMLElement[], options: ElementOptions | null = null): HTMLDivElement {
  let c = block(items, options);
  $(c).addClass('ui-block-25');
  return c;
}

export function block50(items: HTMLElement[], options: ElementOptions | null = null): HTMLDivElement {
  let c = block(items, options);
  $(c).addClass('ui-block-50');
  return c;
}

export function p(text: string, options: any = null): HTMLParagraphElement {
  let c = document.createElement('p');
  c.textContent = text;
  $(c).addClass('ui-p');
  _options(c, options);
  return c;
}

export function panel(items: HTMLElement[] = [], options?: string | ElementOptions): HTMLDivElement {
  let e = div(items, options);
  $(e).addClass('ui-panel');
  return e;
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
    $(d).append(children.map(x => render(x)));
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

export function wideForm(children: InputBase[] = [], options: {} | null = null): HTMLDivElement {
  let d = form(children, options);
  $(d).addClass('ui-form-wide');
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
    $(e).append(children.map(x => render(x)));
  l.textContent = ' ';
  $(l).addClass('ui-label ui-input-label');
  $(d).addClass('ui-input-root ui-input-buttons').append(l).append(e);
  return d;
}

export function contextActions(x: any): HTMLElement {
  return api.grok_UI_Context_Actions(x);
}

export function star(id: string): HTMLElement {
  return api.grok_UI_Star(id);
}

function _icon(type: string, handler: Function, tooltipMsg: string | null = null): HTMLElement {
  let e = $(`<i class="grok-icon grok-font-icon-${type}"></i>`)[0] as HTMLElement;
  e?.addEventListener('click', (e) => handler(e));
  tooltip.bind(e, tooltipMsg);
  return e;
}

function _iconFA(type: string, handler: Function | null, tooltipMsg: string | null = null): HTMLElement {
  let e = $(`<i class="grok-icon fal fa-${type}"></i>`)[0] as HTMLElement;
  if (handler != null)
    e?.addEventListener('click', (e) => handler(e));
  tooltip.bind(e, tooltipMsg);
  return e;
}

export let icons = {
  close: (handler: Function, tooltipMsg: string | null = null) => _icon('close', handler, tooltipMsg),
  help: (handler: Function, tooltipMsg: string | null = null) => _icon('help', handler, tooltipMsg),
  settings: (handler: Function, tooltipMsg: string | null = null) => _icon('settings', handler, tooltipMsg),
  edit: (handler: Function, tooltipMsg: string | null = null) => _iconFA('pen', handler, tooltipMsg),
  save: (handler: Function | null, tooltipMsg: string | null = null) => _iconFA('save', handler, tooltipMsg),
  copy: (handler: Function, tooltipMsg: string | null = null) => _iconFA('copy', handler, tooltipMsg),
  add: (handler: Function, tooltipMsg: string | null = null) => _iconFA('plus', handler, tooltipMsg),
  remove: (handler: Function, tooltipMsg: string | null = null) => _iconFA('minus', handler, tooltipMsg),
  delete: (handler: Function, tooltipMsg: string | null = null) => _iconFA('trash-alt', handler, tooltipMsg),
  undo: (handler: Function, tooltipMsg: string | null = null) => _iconFA('undo', handler, tooltipMsg),
  sync: (handler: Function, tooltipMsg: string | null = null) => _iconFA('sync', handler, tooltipMsg),
  info: (handler: Function, tooltipMsg: string | null = null) => _iconFA('info-circle', handler, tooltipMsg),
  search: (handler: Function, tooltipMsg: string | null = null) => _iconFA('search', handler, tooltipMsg),
  filter: (handler: Function, tooltipMsg: string | null = null) => _iconFA('filter', handler, tooltipMsg),
  play: (handler: Function, tooltipMsg: string | null = null) => _iconFA('play', handler, tooltipMsg),
}

export function setDisplayAll(elements: HTMLElement[], show: boolean): void {
  elements.forEach((e) => setDisplay(e, show));
}

export function setDisplay(element: HTMLElement, show: boolean) {
  if (show)
    element.style.removeProperty('display');
  else
    element.style.display = 'none';
  return element;
}

export function fileBrowser(params: {path?: string, dataSourceFilter?: fileShares[]} = {}): Widget {
  return FilesWidget.create(params);
}

export namespace tools {
  export function click<T extends HTMLElement>(e: T, handler: () => void): T {
    e.addEventListener('click', handler);
    return e;
  }
}

export namespace cards {

  /** Two columns, with picture on the left and details on the right */
  export function summary(picture: HTMLElement, details: any[]): HTMLElement {
    return divH([
      $(picture).addClass('ui-card-picture').get()[0],
      div(details.map((d) => render(d)), 'ui-card-details')
    ], 'ui-card')
  }
}

export namespace labels {

  /** Minor details to show using the smaller font */
  export function details(s: string): HTMLDivElement {
    return divText(s, 'ui-label-details');
  }
}

/** Visual hints attached to an element. They can be used to introduce a new app to users. */
export namespace hints {

  /** Adds a hint indication to the provided element and returns it. */
  export function addHint(el: HTMLElement): HTMLElement {
    el.classList.add('ui-hint-target');
    const hintIndicator = div([], 'ui-hint-blob');
    el.append(hintIndicator);

    const width = el ? el.clientWidth : 0;
    const height = el ? el.clientHeight : 0;

    const hintNode = el.getBoundingClientRect();
    const indicatorNode = hintIndicator.getBoundingClientRect();

    const hintPosition = $(el).css('position') ?? 'static';

    function setDefaultStyles() {
      $(hintIndicator).css('position', 'absolute');
      $(hintIndicator).css('left', '0');
      $(hintIndicator).css('top', '0');
    };

    const positions = {
      static: () => {
        $(hintIndicator).css('position', 'absolute');
        $(hintIndicator).css('margin-left', hintNode.left + 1 == indicatorNode.left ? 0 : -width);
        $(hintIndicator).css('margin-top', hintNode.top + 1 == indicatorNode.top ? 0 : -height);
      },
      relative: setDefaultStyles,
      fixed: setDefaultStyles,
      absolute: setDefaultStyles,
      sticky: setDefaultStyles,
    };

    positions[hintPosition]();

    if ($(el).hasClass('d4-ribbon-item')){
      $(hintIndicator).css('margin-left', '0px')
    }

    return el;
  }

  /** Removes the hint indication from the provided element and returns it. */
  export function remove(el: HTMLElement): HTMLElement {
    $(el).find('div.ui-hint-blob')[0]?.remove();
    el.classList.remove('ui-hint-target');
    return el;
  }
}

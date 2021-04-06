/**
 * Routines for building UI
 * @module ui
 **/
var __createBinding = (this && this.__createBinding) || (Object.create ? (function(o, m, k, k2) {
  if (k2 === undefined) k2 = k;
  Object.defineProperty(o, k2, { enumerable: true, get: function() { return m[k]; } });
}) : (function(o, m, k, k2) {
  if (k2 === undefined) k2 = k;
  o[k2] = m[k];
}));
var __setModuleDefault = (this && this.__setModuleDefault) || (Object.create ? (function(o, v) {
  Object.defineProperty(o, 'default', { enumerable: true, value: v });
}) : function(o, v) {
  o['default'] = v;
});
var __importStar = (this && this.__importStar) || function (mod) {
  if (mod && mod.__esModule) return mod;
  var result = {};
  if (mod != null) for (var k in mod) if (k !== 'default' && Object.prototype.hasOwnProperty.call(mod, k)) __createBinding(result, mod, k);
  __setModuleDefault(result, mod);
  return result;
};
var __importDefault = (this && this.__importDefault) || function (mod) {
  return (mod && mod.__esModule) ? mod : { 'default': mod };
};
define(['require', 'exports', './src/view', './src/widgets', './src/wrappers', './src/functions', 'cash-dom', './src/events', './src/utils', 'rxjs', './src/grid'], function (require, exports, view_1, widgets_1, wrappers_1, functions_1, cash_dom_1, events_1, utils_1, rxjs, grid_1) {
  'use strict';
  Object.defineProperty(exports, '__esModule', { value: true });
  exports.buttonsInput = exports.narrowForm = exports.form = exports.label = exports.panel = exports.p = exports.block50 = exports.block25 = exports.block75 = exports.block = exports.splitH = exports.splitV = exports.boxFixed = exports.box = exports.EntityMetaDartProxy = exports.ObjectHandler = exports.ObjectHandlerResolutionArgs = exports.tooltip = exports.Tooltip = exports.tools = exports.onSizeChanged = exports.textInput = exports.columnsInput = exports.columnInput = exports.moleculeInput = exports.boolInput = exports.dateInput = exports.floatInput = exports.stringInput = exports.multiChoiceInput = exports.choiceInput = exports.intInput = exports.tree = exports.inputs = exports.popupMenu = exports.virtualView = exports.rangeSlider = exports.showPopup = exports.bind = exports.dialog = exports.list = exports.wait = exports.table = exports.tableFromMap = exports.comboPopupItems = exports.comboPopup = exports.bigButton = exports.button = exports.setUpdateIndicator = exports.loader = exports.card = exports.divH = exports.divV = exports.div = exports.inlineText = exports.renderInline = exports.span = exports.renderCard = exports.render = exports.extract = exports.iconFA = exports.divText = exports.tabControl = exports.accordion = exports.h3 = exports.h2 = exports.h1 = exports.canvas = exports._backColor = exports._color = exports._class = exports._innerText = exports.empty = exports.appendAll = exports.element = void 0;
  cash_dom_1 = __importDefault(cash_dom_1);
  rxjs = __importStar(rxjs);
  let api = window;
  /**
     * Creates an instance of the element for the specified tag, and optionally assigns it a CSS class.
     * @param {string} tagName The name of an element.
     * @param {string | null} className
     * @returns {HTMLElement}
     */
  function element(tagName, className = null) {
    let x = document.createElement(tagName);
    if (className !== null)
      _class(x, className);
    return x;
  }
  exports.element = element;
  /** Appends multiple elements to root, and returns root.
     *  An element could be either {@link HTMLElement} or {@link Viewer}.
     *
     * @param {HTMLElement} root
     * @param {(HTMLElement | Viewer)[]} elements
     * @returns {HTMLElement} */
  function appendAll(root, elements) {
    let fragment = document.createDocumentFragment();
    for (let e of elements)
      fragment.appendChild(render(e));
    root.appendChild(fragment);
    return root;
  }
  exports.appendAll = appendAll;
  function empty(e) {
    while (e.firstChild)
      e.removeChild(e.firstChild);
    return e;
  }
  exports.empty = empty;
  function _innerText(x, s) {
    x.innerText = s;
    return x;
  }
  exports._innerText = _innerText;
  function _class(x, s) {
    x.classList.add(s);
    return x;
  }
  exports._class = _class;
  function _color(x, s) {
    x.style.color = s;
    return x;
  }
  exports._color = _color;
  function _backColor(x, s) {
    x.style.backgroundColor = s;
    return x;
  }
  exports._backColor = _backColor;
  /**
     * @param {number} height
     * @param {number} width
     * @returns {HTMLCanvasElement}
     * */
  function canvas(width = null, height = null) {
    let result = element('CANVAS');
    if (height != null && width != null) {
      cash_dom_1.default(result).height(height);
      cash_dom_1.default(result).width(width);
      cash_dom_1.default(result).prop('height', `${height}`);
      cash_dom_1.default(result).prop('width', `${width}`);
    }
    return result;
  }
  exports.canvas = canvas;
  /** @returns {HTMLHeadingElement} */
  function h1(s) {
    return _innerText(element('h1'), s);
  }
  exports.h1 = h1;
  /** @returns {HTMLHeadingElement} */
  function h2(s) {
    let x = element('h2');
    x.innerText = s;
    return x;
  }
  exports.h2 = h2;
  /** @returns {HTMLHeadingElement} */
  function h3(s) {
    let x = element('h3');
    x.innerText = s;
    return x;
  }
  exports.h3 = h3;
  /** @returns {Accordion} */
  function accordion(key = null) {
    return widgets_1.Accordion.create(key);
  }
  exports.accordion = accordion;
  /**
     * @param {Object} pages - list of page factories
     * @param {boolean} vertical
     * @returns {TabControl} */
  function tabControl(pages = null, vertical = false) {
    let tabs = widgets_1.TabControl.create(vertical);
    if (pages != null) {
      for (let key of Object.keys(pages)) {
        let value = pages[key];
        tabs.addPane(key, value instanceof Function ? value : () => render(value));
      }
    }
    return tabs;
  }
  exports.tabControl = tabControl;
  /** Returns DivElement with the specified inner text
     * @param {string} text
     * @param {string | ElementOptions | null} options
     * @returns {HTMLDivElement} */
  function divText(text, options = null) {
    let e = element('div');
    e.innerText = text;
    utils_1._options(e, options);
    return e;
  }
  exports.divText = divText;
  /** Returns a font-awesome icon with the specified name, handler, and tooltip.
     * Sample: {@link https://public.datagrok.ai/js/samples/ui/icons}
     * @param {string} name - icon name (omit the "fa-" prefix)
     * @param {Function} handler
     * @param {String} tooltipMsg
     * @returns {HTMLElement} */
  function iconFA(name, handler, tooltipMsg = null) {
    let i = element('i');
    i.classList.add('grok-icon');
    i.classList.add('fal');
    i.classList.add(`fa-${name}`);
    i.addEventListener('click', handler);
    if (tooltipMsg !== null)
      exports.tooltip.bind(i, tooltipMsg);
    return i;
  }
  exports.iconFA = iconFA;
  function extract(x) {
    if (x == null)
      return null;
    if (x instanceof widgets_1.Widget)
      return x.root;
    if (typeof x.root !== 'undefined')
      return x.root;
    if (typeof x.wrap === 'function') {
      if (x.length === 1)
        return x[0];
      else
        return span(x.get());
    } // jquery (cash) object
    return x;
  }
  exports.extract = extract;
  /** Renders object to html element.
     * @param {object} x
     * @returns {HTMLElement} */
  function render(x) {
    x = extract(x);
    // @ts-ignore
    if (api.grok_UI_Render == null)
      return x;
    return api.grok_UI_Render(x);
  }
  exports.render = render;
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
  function renderCard(x) {
    return api.grok_UI_RenderCard(x);
  }
  exports.renderCard = renderCard;
  /** Renders span
     * @param {object[]} x
     * @returns {HTMLElement}. */
  function span(x) {
    return api.grok_UI_Span(x);
  }
  exports.span = span;
  function renderInline(x) {
    let handler = ObjectHandler.forEntity(x);
    return handler == null ? render(x) : handler.renderMarkup(x);
  }
  exports.renderInline = renderInline;
  /** Renders inline text, calling [renderMarkup] for each non-HTMLElement
     * @param {object[]} objects
     * @returns {HTMLElement}. */
  function inlineText(objects) {
    return span(objects.map((item) => renderInline(item)));
  }
  exports.inlineText = inlineText;
  /**
     * @param {object[]} children
     * @param {string | ElementOptions} options
     * @returns {HTMLDivElement}
     * */
  function div(children = [], options = null) {
    if (!Array.isArray(children))
      children = [children];
    let d = document.createElement('div');
    if (children != null) {
      cash_dom_1.default(d).append(children.map(render));
    }
    utils_1._options(d, options);
    cash_dom_1.default(d).addClass('ui-div');
    return d;
  }
  exports.div = div;
  /** Div flex-box container that positions child elements vertically.
     * @param {object[]} items
     * @param {string | ElementOptions} options
     * @returns {HTMLDivElement} */
  function divV(items, options = null) {
    return utils_1._options(api.grok_UI_DivV(items == null ? null : items.map(render), null), options);
  }
  exports.divV = divV;
  /** Div flex-box container that positions child elements horizontally.
     * @param {object[]} items
     * @param {string | ElementOptions} options
     * @returns {HTMLDivElement} */
  function divH(items, options = null) {
    return utils_1._options(api.grok_UI_DivH(items == null ? null : items.map(render), null), options);
  }
  exports.divH = divH;
  /** Renders content as a card. */
  function card(content) {
    return div([content], 'd4-item-card');
  }
  exports.card = card;
  function loader() {
    return api.grok_UI_Loader();
  }
  exports.loader = loader;
  function setUpdateIndicator(element, updating = true) {
    return api.grok_UI_SetUpdateIndicator(element, updating);
  }
  exports.setUpdateIndicator = setUpdateIndicator;
  /**
     * Creates a button with the specified text, click handler, and tooltip
     * @param {string | Element | Array<string | Element>} content
     * @param {Function} handler
     * @param {string} tooltip
     * @returns {HTMLButtonElement}
     * */
  function button(content, handler, tooltip = null) {
    return api.grok_UI_Button(content, handler, tooltip);
  }
  exports.button = button;
  function bigButton(text, handler, tooltip = null) {
    return api.grok_UI_BigButton(text, handler, tooltip);
  }
  exports.bigButton = bigButton;
  /**
     * Creates a combo popup with the specified icons and items
     * @param {string | HTMLElement} caption
     * @param {Array<string>} items
     * @param {Function} handler (item) => {...}
     * @param {Function} renderer (item) => {...}
     * @returns {HTMLElement}
     * */
  function comboPopup(caption, items, handler, renderer = null) {
    return api.grok_UI_ComboPopup(caption, items, handler, renderer !== null ? (item) => renderer(wrappers_1.toJs(item)) : null);
  }
  exports.comboPopup = comboPopup;
  /**
     * Creates a combo popup with the specified icons and items
     * @param {string | HTMLElement} caption
     * @param {Map<string>} items
     * @returns {HTMLElement}
     * */
  function comboPopupItems(caption, items) {
    return api.grok_UI_ComboPopup(caption, Object.keys(items), (key) => items[key]());
  }
  exports.comboPopupItems = comboPopupItems;
  /** Creates a visual table based on [map]. */
  function tableFromMap(map) {
    return api.grok_UI_TableFromMap(map);
  }
  exports.tableFromMap = tableFromMap;
  /** Creates a visual table based on [items] and [renderer]. */
  function table(items, renderer, columnNames = null) {
    return api.grok_HtmlTable(items, renderer !== null ? (object, ind) => renderer(wrappers_1.toJs(object), ind) : null, columnNames).root;
  }
  exports.table = table;
  /** Waits for Future<Element> function to complete and collect its result.*/
  function wait(getElement) {
    return wrappers_1.toJs(api.grok_UI_Wait(getElement));
  }
  exports.wait = wait;
  /** Creates a visual element representing list of [items]. */
  function list(items) {
    return api.grok_UI_List(Array.from(items).map(wrappers_1.toDart));
  }
  exports.list = list;
  /** Creates a [Dialog].
     * @returns {Dialog} */
  function dialog(title = '') {
    return widgets_1.Dialog.create(title);
  }
  exports.dialog = dialog;
  /** Binds [item] with the [element]. It enables selecting it as a current object, drag-and-drop,
     * tooltip, and popup menu.
     * @param item
     * @param {Element} element
     * @returns {Element}. */
  function bind(item, element) {
    return api.grok_UI_Bind(item, element);
  }
  exports.bind = bind;
  /** Shows popup with the [element] near the [anchor].
     * tooltip, and popup menu.
     * @param {Element} element
     * @param {Element} anchor
     * @param {boolean} vertical
     * @returns {Element}. */
  function showPopup(element, anchor, vertical = false) {
    return api.grok_UI_ShowPopup(element, anchor, vertical);
  }
  exports.showPopup = showPopup;
  function rangeSlider(minRange, maxRange, min, max) {
    let rs = widgets_1.RangeSlider.create();
    rs.setValues(minRange, maxRange, min, max);
    return rs;
  }
  exports.rangeSlider = rangeSlider;
  /**
     * Creates a virtual list widget
     * @param {number} length - number of elements
     * @param {Function} renderer
     * @param {boolean} verticalScroll - vertical or horizontal scrolling
     * @param {number} maxColumns - maximum number of items on the non-scrolling axis
     * @returns {VirtualView}
     */
  function virtualView(length, renderer, verticalScroll = true, maxColumns = 1000) {
    let view = view_1.VirtualView.create(verticalScroll, maxColumns);
    view.setData(length, renderer);
    return view;
  }
  exports.virtualView = virtualView;
  function popupMenu(items) {
    function populate(menu, item) {
      for (let key of Object.keys(item)) {
        let value = item[key];
        if (value instanceof Function)
          menu.item(key, value);
        else
          populate(menu.group(key), value);
      }
    }
    let menu = widgets_1.Menu.popup();
    populate(menu, items);
    menu.show();
  }
  exports.popupMenu = popupMenu;
  function inputs(inputs, options) {
    return form(inputs, options);
  }
  exports.inputs = inputs;
  /** Creates new nodes tree
     * @returns {TreeViewNode} */
  function tree() {
    return widgets_1.TreeViewNode.tree();
  }
  exports.tree = tree;
  function intInput(name, value, onValueChanged = null) {
    return new widgets_1.InputBase(api.grok_IntInput(name, value), onValueChanged);
  }
  exports.intInput = intInput;
  function choiceInput(name, selected, items, onValueChanged = null) {
    return new widgets_1.InputBase(api.grok_ChoiceInput(name, selected, items), onValueChanged);
  }
  exports.choiceInput = choiceInput;
  function multiChoiceInput(name, value, items, onValueChanged = null) {
    return new widgets_1.InputBase(api.grok_MultiChoiceInput(name, value, items), onValueChanged);
  }
  exports.multiChoiceInput = multiChoiceInput;
  function stringInput(name, value, onValueChanged = null) {
    return new widgets_1.InputBase(api.grok_StringInput(name, value), onValueChanged);
  }
  exports.stringInput = stringInput;
  function floatInput(name, value, onValueChanged = null) {
    return new widgets_1.InputBase(api.grok_FloatInput(name, value), onValueChanged);
  }
  exports.floatInput = floatInput;
  function dateInput(name, value, onValueChanged = null) {
    return new widgets_1.InputBase(api.grok_DateInput(name, value.d), onValueChanged);
  }
  exports.dateInput = dateInput;
  function boolInput(name, value, onValueChanged = null) {
    return new widgets_1.InputBase(api.grok_BoolInput(name, value), onValueChanged);
  }
  exports.boolInput = boolInput;
  function moleculeInput(name, value, onValueChanged = null) {
    return new widgets_1.InputBase(api.grok_MoleculeInput(name, value), onValueChanged);
  }
  exports.moleculeInput = moleculeInput;
  function columnInput(name, table, value, onValueChanged = null) {
    return new widgets_1.InputBase(api.grok_ColumnInput(name, table.d, value.d), onValueChanged);
  }
  exports.columnInput = columnInput;
  function columnsInput(name, table, onValueChanged = null) {
    return new widgets_1.InputBase(api.grok_ColumnsInput(name, table.d), onValueChanged);
  }
  exports.columnsInput = columnsInput;
  function textInput(name, value, onValueChanged = null) {
    return new widgets_1.InputBase(api.grok_TextInput(name, value), onValueChanged);
  }
  exports.textInput = textInput;
  /**
     * @param {HTMLElement} element
     * @returns {rxjs.Observable} */
  function onSizeChanged(element) {
    if (utils_1._isDartium()) {
      // Polyfill for Dartium which does not have a native ResizeObserver
      return new rxjs.Observable(function (observer) {
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
    return rxjs.Observable.create(function (observer) {
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
  exports.onSizeChanged = onSizeChanged;
  /** UI Tools **/
  class tools {
    static handleResize(element, onChanged) {
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
  exports.tools = tools;
  /** Represents a tooltip. */
  class Tooltip {
    /** Hides the tooltip. */
    hide() {
      api.grok_Tooltip_Hide();
    }
    /** Associated the specified visual element with the corresponding item. */
    bind(element, tooltip) {
      api.grok_Tooltip_SetOn(element, tooltip);
      return element;
    }
    /** Shows the tooltip at the specified position
         * @param {HTMLElement | string} content
         * @param {number} x
         * @param {number} y */
    show(content, x, y) {
      api.grok_Tooltip_Show(content, x, y);
    }
    showRowGroup(dataFrame, indexPredicate, x, y) {
      api.grok_Tooltip_ShowRowGroup(dataFrame.d, indexPredicate, x, y);
    }
    /** Returns a tooltip element.
         * @returns {HTMLElement} */
    get root() {
      return api.grok_Tooltip_Get_Root();
    }
    /** @returns {boolean} isVisible */
    get isVisible() {
      return api.grok_Tooltip_Is_Visible();
    }
    /** @returns {Observable} */
    get onTooltipRequest() {
      return events_1.__obs('d4-tooltip-request');
    }
    /** @returns {Observable} */
    get onTooltipShown() {
      return events_1.__obs('d4-tooltip-shown');
    }
    /** @returns {Observable} */
    get onTooltipClosed() {
      return events_1.__obs('d4-tooltip-closed');
    }
  }
  exports.Tooltip = Tooltip;
  exports.tooltip = new Tooltip();
  class ObjectHandlerResolutionArgs {
    constructor(object, context) {
      this.object = object;
      this.context = context;
      this.handler = null;
    }
  }
  exports.ObjectHandlerResolutionArgs = ObjectHandlerResolutionArgs;
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
  class ObjectHandler {
    /** Type of the object that this meta handles. */
    get type() {
      throw 'Not defined.';
    }
    /**
         * Override this method to check whether this meta class should handle the specified object.
         * @param x - specified object.
         * @returns {boolean}
         * */
    isApplicable(x) {
      throw 'Not defined.';
    }
    /** String representation of the [item], by default item.toString().
         * @param x - item
         * @returns {string} */
    getCaption(x) {
      return `${x}`;
    }
    /** @returns {CanvasRenderer} */
    getCanvasRenderer() {
      return null;
    }
    /** @returns {GridCellRenderer} */
    getGridCellRenderer() {
      return null;
    }
    /** Renders icon for the item.
         * @param x - item
         * @returns {Element} */
    renderIcon(x, context = null) {
      return divText(this.getCaption(x));
    }
    /** Renders markup for the item.
         * @param x - item
         * @returns {Element} */
    renderMarkup(x, context = null) {
      return divText(this.getCaption(x));
    }
    /** Renders tooltip for the item.
         * @param x - item
         * @returns {Element} */
    renderTooltip(x, context = null) {
      return divText(this.getCaption(x));
    }
    /** Renders card div for the item.
         * @param x - item
         * @returns {Element} */
    renderCard(x, context = null) {
      return divText(this.getCaption(x));
    }
    /** Renders properties list for the item.
         * @param x - item
         * @returns {Element} */
    renderProperties(x, context = null) {
      return divText(this.getCaption(x));
    }
    /** Renders view for the item.
         * @param x - item
         * @returns {Element} */
    renderView(x, context = null) {
      return this.renderProperties(x);
    }
    /** Gets called once upon the registration of meta export class. */
    init() { }
    /** Registers entity handler.
         * @param {ObjectHandler} meta */
    static register(meta) {
      api.grok_Meta_Register(meta);
      let cellRenderer = meta.getGridCellRenderer();
      if (cellRenderer != null)
        grid_1.GridCellRenderer.register(cellRenderer);
    }
    /** @returns {ObjectHandler[]} */
    static list() {
      return wrappers_1.toJs(api.grok_Meta_List());
    }
    static onResolve(observer) {
      return _objectHandlerSubject.subscribe(observer);
    }
    /**
         * @param {Object} object
         * @param {Object} context
         * @returns {ObjectHandler}
         * */
    static forEntity(object, context = null) {
      let args = new ObjectHandlerResolutionArgs(object, context);
      _objectHandlerSubject.next(args);
      if (args.handler !== null)
        return args.handler;
      return wrappers_1.toJs(api.grok_Meta_ForEntity(wrappers_1.toDart(object)));
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
    registerParamFunc(name, run) {
      // @ts-ignore
      new functions_1.Functions().registerParamFunc(name, this.type, run, this.isApplicable);
    }
  }
  exports.ObjectHandler = ObjectHandler;
  class EntityMetaDartProxy extends ObjectHandler {
    constructor(d) {
      super();
      this.d = d;
    }
    get type() { return api.grok_Meta_Get_Type(this.d); }
    isApplicable(x) { return api.grok_Meta_IsApplicable(this.d, wrappers_1.toDart(x)); }
    getCaption(x) { return api.grok_Meta_Get_Name(this.d, x); }
    renderIcon(x, context = null) { return api.grok_Meta_RenderIcon(x); }
    renderMarkup(x, context = null) { return api.grok_Meta_RenderMarkup(x); }
    renderTooltip(x, context = null) { return api.grok_Meta_RenderTooltip(x); }
    renderCard(x, context = null) { return api.grok_Meta_RenderCard(x); }
    renderProperties(x, context = null) { return api.grok_Meta_RenderProperties(x); }
    renderView(x, context = null) { return api.grok_Meta_RenderProperties(x); }
  }
  exports.EntityMetaDartProxy = EntityMetaDartProxy;
  function box(item, options = null) {
    if (item instanceof widgets_1.Widget) {
      item = item.root;
    }
    if (item instanceof widgets_1.InputBase) {
      console.log('inputbase');
      item = item.root;
    }
    if (item != null && cash_dom_1.default(item).hasClass('ui-box')) {
      return item;
    }
    let c = document.createElement('div');
    utils_1._options(c, options);
    cash_dom_1.default(c).append(item).addClass('ui-box');
    return c;
  }
  exports.box = box;
  function boxFixed(item, options = null) {
    let c = box(item, options);
    cash_dom_1.default(c).addClass('ui-box-fixed');
    return c;
  }
  exports.boxFixed = boxFixed;
  /** Div flex-box container that positions child elements vertically.
     *  @param {object[]} items */
  function splitV(items, options = null) {
    let b = box(null, options);
    cash_dom_1.default(b).addClass('ui-split-v').append(items.map(item => box(item)));
    return b;
  }
  exports.splitV = splitV;
  /** Div flex-box container that positions child elements horizontally.
     *  @param {object[]} items */
  function splitH(items, options = null) {
    let b = box(null, options);
    cash_dom_1.default(b).addClass('ui-split-h').append(items.map(item => box(item)));
    return b;
  }
  exports.splitH = splitH;
  function block(items, options = null) {
    let c = div(items, options);
    cash_dom_1.default(c).addClass('ui-block');
    return c;
  }
  exports.block = block;
  function block75(items, options = null) {
    return cash_dom_1.default(block(items, options)).addClass('ui-block-75')[0];
  }
  exports.block75 = block75;
  function block25(items, options = null) {
    return cash_dom_1.default(block(items, options)).addClass('ui-block-25')[0];
  }
  exports.block25 = block25;
  function block50(items, options = null) {
    return cash_dom_1.default(block(items, options)).addClass('ui-block-50')[0];
  }
  exports.block50 = block50;
  function p(text, options = null) {
    let c = document.createElement('p');
    c.textContent = text;
    cash_dom_1.default(c).addClass('ui-p');
    utils_1._options(c, options);
    return c;
  }
  exports.p = p;
  function panel(items, options) {
    return cash_dom_1.default(div(items, options)).addClass('ui-panel')[0];
  }
  exports.panel = panel;
  function label(text, options = null) {
    let c = document.createElement('label');
    c.textContent = text;
    cash_dom_1.default(c).addClass('ui-label');
    utils_1._options(c, options);
    return c;
  }
  exports.label = label;
  function form(children = [], options = null) {
    let d = document.createElement('div');
    if (children != null) {
      cash_dom_1.default(d).append(children.map(render));
    }
    utils_1._options(d, options);
    cash_dom_1.default(d).addClass('ui-form');
    return d;
  }
  exports.form = form;
  function narrowForm(children = [], options = null) {
    let d = form(children, options);
    cash_dom_1.default(d).addClass('ui-form-condensed');
    return d;
  }
  exports.narrowForm = narrowForm;
  function buttonsInput(children = []) {
    if (!Array.isArray(children))
      children = [children];
    let d = document.createElement('div');
    let l = document.createElement('label');
    let e = document.createElement('div');
    cash_dom_1.default(e).addClass('ui-input-editor');
    if (children != null)
      cash_dom_1.default(e).append(children.map(render));
    l.textContent = ' ';
    cash_dom_1.default(l).addClass('ui-label ui-input-label');
    cash_dom_1.default(d).addClass('ui-input-root ui-input-buttons').append(l).append(e);
    return d;
  }
  exports.buttonsInput = buttonsInput;
});

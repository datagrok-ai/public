/**
 * Routines for building UI
 * @module ui
 **/

import {Viewer} from "./src/viewer";
import {VirtualView} from "./src/view";
import {Accordion, Dialog, InputBase, Menu, TabControl, TreeViewNode, Widget} from "./src/widgets";
import {toDart} from "./src/wrappers";

/**
 * Creates an instance of the element for the specified tag, and optionally assigns it a CSS class.
 * @param {string} tagName The name of an element.
 * @param {string | null} className
 */
export function element(tagName, className = null) {
        let x = document.createElement(tagName);
        if (className !== null)
            _class(x, className);
        return x;
    }

/** Appends multiple elements to root, and returns root.
 *  An element could be either {@link HTMLElement} or {@link Viewer}.
 *
 * @param {HTMLElement} root
 * @param {object[]} elements
 * @returns {HTMLElement} */
export function appendAll(root, elements) {
        let fragment = document.createDocumentFragment();
        for (let i = 0; i < elements.length; i++) {
            let e = elements[i];
            if (e instanceof Viewer)
                e = e.root;
            fragment.appendChild(e);
        }
        root.appendChild(fragment);
        return root;
    }

export function empty(e) {
    while (e.firstChild)
        e.removeChild(foo.firstChild);
    return e;
}

export function _innerText(x, s) { x.innerText = s; return x; }
export function _class(x, s) { x.classList.add(s); return x; }
export function _color(x, s) { x.style.color = s; return x; }
export function _backColor(x, s) { x.style.backgroundColor = s; return x; }

/** @returns {HTMLCanvasElement} */
export function canvas() { return element("CANVAS"); }

/** @returns {HTMLHeadingElement} */
export function h1(s) { return _innerText(ui.element('h1'), s); }

/** @returns {HTMLHeadingElement} */
export function h2(s) { let x = element('h2'); x.innerText = s; return x; }

/** @returns {HTMLHeadingElement} */
export function h3(s) { let x = element('h3'); x.innerText = s; return x; }

/** @returns {Accordion} */
export function accordion() { return Accordion.create(); }

/** @returns {TabControl}
 * @param {boolean} vertical */
export function tabControl(vertical = false) { return TabControl.create(vertical); }

/** Returns DivElement with the specified inner text
 * @param {string} text
 * @returns {HTMLDivElement} */
export function divText(text, className = null) {
    let e = element('div');
    e.innerText = text;
    if (className != null)
        e.classList.add(className);

    return e;
}

/** Returns a font-awesome icon with the specified name, handler, and tooltip.
 * Sample: {@link https://public.datagrok.ai/js/samples/ui/icons}
 * @param {string} name - icon name (omit the "fa-" prefix)
 * @param {Function} handler
 * @param {String} tooltipMsg
 * @returns {HTMLElement} */
export function iconFA(name, handler, tooltipMsg = null) {
    let i = element('i');
    i.classList.add('grok-icon');
    i.classList.add('fal');
    i.classList.add(`fa-${name}`);
    i.addEventListener("click", handler);
    if (tooltipMsg !== null)
        tooltip(i, tooltipMsg);
    return i;
}

/** Renders object to html element.
 * @param {object} x
 * @returns {HTMLElement} */
export function render(x) {
    if (x instanceof Widget)
        return x.root;
    return grok_UI_Render(x);
}

/** Renders a table cell to html element, taking into account grid's formatting and color-coding.
 * @param {Cell} tableCell
 * @returns {HTMLElement} */
export function renderCell(tableCell) {

}

/** Renders a table cell to html element,
 * Takes into account grid's formatting and color-coding.
 * @param {Row} table
 * @param {string[]} columnNames
 * @returns {HTMLElement} */
export function renderForm(table, columnNames = null) {
  return null;
}

/** Renders object to html card.
 * @param {object} x
 * @returns {HTMLElement}. */
export function renderCard(x) { return grok_UI_RenderCard(x); }

/** Renders span
 * @param {object[]} x
 * @returns {HTMLElement}. */
export function span(x) { return grok_UI_Span(x); }

/** @returns {HTMLDivElement} */
export function div(items = [], className = null) { return grok_UI_Div(items, className); }

/** Div flex-box container that positions child elements vertically.
 *  @param {HTMLElement[]} items
 *  @param {string} className - comma-separated CSS class names
 *  @returns {HTMLDivElement} */
export function divV(items, className = null) { return grok_UI_DivV(items, className); }

/** Div flex-box container that positions child elements horizontally.
 *  @param {HTMLElement[]} items
 *  @param {string} className - comma-separated CSS class names
 *  @returns {HTMLDivElement} */
export function divH(items, className = null) { return grok_UI_DivH(items, className); }

/** Renders content as a card. */
export function card(content) { return ui.div([content], 'd4-item-card'); }

export function loader() { return grok_UI_Loader(); }
export function setUpdateIndicator(element, updating = true) { return grok_UI_SetUpdateIndicator(element, updating)};

/**
 * Creates a button with the specified text, click handler, and tooltip
 * @param {string} text
 * @param {Function} handler
 * @param {string} tooltip
 * @returns {HTMLButtonElement}
 * */
export function button(text, handler, tooltip = null) { return grok_UI_Button(text, handler, tooltip); }

export function bigButton(text, handler, tooltip = null) { return grok_UI_BigButton(text, handler, tooltip); }

/**
 * Creates a combo popup with the specified icons and items
 * @param {HTMLElement} icon
 * @param {Array<string>} items
 * @param {Function} handler (item) => {...}
 * @returns {HTMLElement}
 * */
export function comboPopup(icon, items, handler) { return grok_UI_ComboPopup(icon, items, handler); }

/** Creates a visual table based on [map]. */
export function tableFromMap(map) { return grok_UI_TableFromMap(map); }

/** Creates a visual element representing list of [items]. */
export function list(items) { return grok_UI_List(Array.from(items).map(toDart)); }

/** Creates a [Dialog].
 * @returns {Dialog} */
export function dialog(title = '') { return Dialog.create(title); }

/** Binds [item] with the [element]. It enables selecting it as a current object, drag-and-drop,
 * tooltip, and popup menu.
 * @returns Element. */
export function bind(item, element) { return grok_UI_Bind(item, element); }

/**
 * Creates a virtual list widget
 * @param {number} length - number of elements
 * @param {Function} renderer
 * @returns {VirtualView}
 */
export function virtualView(length, renderer) {
    let view = VirtualView.create();
    view.setData(length, renderer);
    return view;
}

export function popupMenu(items) {
    function populate(menu, item) {
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

export function tooltipHide() { grok_Tooltip_Hide(); }

export function tooltip(e, x) { grok_Tooltip_SetOn(e, x); return e; }

export function tooltipShow(content, x, y) { grok_Tooltip_Show(content, x, y); }

export function inputs(inputs) { return div(inputs.map((x) => x.root), 'pure-form,pure-form-aligned');}

/** Creates new nodes tree
 * @returns {TreeViewNode} */
export function tree() { return TreeViewNode.tree(); }

export function intInput(name, value) { return new InputBase(grok_IntInput(name, value)); }
export function choiceInput(name, selected, items) { return new InputBase(grok_ChoiceInput(name, selected, items)); }
export function multiChoiceInput(name, value, items) { return new InputBase(grok_MultiChoiceInput(name, value, items)); }
export function stringInput(name, value) { return new InputBase(grok_StringInput(name, value)); }
export function floatInput(name, value) { return new InputBase(grok_FloatInput(name, value)); }
export function dateInput(name, value) { return new InputBase(grok_DateInput(name, value.d)); }
export function boolInput(name, value, callback = null) { return new InputBase(grok_BoolInput(name, value), callback); }
export function moleculeInput(name, value) { return new InputBase(grok_MoleculeInput(name, value)); }
export function columnInput(name, table, value) { return new InputBase(grok_ColumnInput(name, table.d, value.d)); }
export function columnsInput(name, table) { return new InputBase(grok_ColumnsInput(name, table.d)); }

/** UI Tools **/
export class tools {

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

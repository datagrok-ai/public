/**
 * Routines for building UI
 * @module ui
 **/

import {Viewer} from "./src/viewer";
import {Accordion, Dialog, InputBase, TabControl, TreeViewNode, Widget} from "./src/widgets";
import {Cell, Row} from "./src/dataframe";
import {VirtualView} from "./src/view";

/**
 * Creates an instance of the element for the specified tag, and optionally assigns it a CSS class.
 * @param {string} tagName The name of an element.
 * @param {string | null} className
 */
export function element(tagName: string, className?: string | null): HTMLElement

/** Appends multiple elements to root, and returns root.
 *  An element could be either {@link HTMLElement} or {@link Viewer}.
 *
 * @param {HTMLElement} root
 * @param {object[]} elements
 * @returns {HTMLElement} */
export function appendAll(root: HTMLElement, elements: (HTMLElement | Viewer)[]): HTMLElement

export function empty<T extends HTMLElement>(e: T): T

/*
export function _innerText(x, s) {
    x.innerText = s;
    return x;
}

export function _class(x, s) {
    x.classList.add(s);
    return x;
}

export function _color(x, s) {
    x.style.color = s;
    return x;
}

export function _backColor(x, s) {
    x.style.backgroundColor = s;
    return x;
}*/

/** @returns {HTMLCanvasElement} */
export function canvas(): HTMLCanvasElement

/** @returns {HTMLHeadingElement} */
export function h1(s: string): HTMLHeadingElement

/** @returns {HTMLHeadingElement} */
export function h2(s: string): HTMLHeadingElement

/** @returns {HTMLHeadingElement} */
export function h3(s: string): HTMLHeadingElement

/** @returns {Accordion} */
export function accordion(): Accordion

/** @returns {TabControl}
 * @param {Object} pages
 * @param {boolean} vertical */
export function tabControl(pages?: any, vertical?: boolean): TabControl

/** Returns DivElement with the specified inner text
 * @param {string} text
 * @returns {HTMLDivElement} */
export function divText(text: string, className?: string | null): HTMLDivElement

/** Returns a font-awesome icon with the specified name, handler, and tooltip.
 * Sample: {@link https://public.datagrok.ai/js/samples/ui/icons}
 * @param {string} name - icon name (omit the "fa-" prefix)
 * @param {Function} handler
 * @param {String} tooltipMsg
 * @returns {HTMLElement} */
export function iconFA(name: string, handler: Function, tooltipMsg?: string | null): HTMLElement

export function extract(x: any) : any

/** Renders object to html element.
 * @param {object} x
 * @returns {HTMLElement} */
export function render(x: Widget | object): HTMLElement

/** Renders a table cell to html element, taking into account grid's formatting and color-coding.
 * @param {Cell} tableCell
 * @returns {HTMLElement} */
export function renderCell(tableCell: Cell): HTMLElement

/** Renders a table cell to html element,
 * Takes into account grid's formatting and color-coding.
 * @param {Row} table
 * @param {string[]} columnNames
 * @returns {HTMLElement} */
export function renderForm(table: Row, columnNames?: string[] | null): HTMLElement | null

/** Renders object to html card.
 * @param {object} x
 * @returns {HTMLElement}. */
export function renderCard(x: object): HTMLElement

/** Renders span
 * @param {object[]} x
 * @returns {HTMLElement}. */
export function span(x: object[]): HTMLElement

/** @returns {HTMLDivElement} */
export function div(items?: any[], className?: string | null): HTMLDivElement

/** Div flex-box container that positions child elements vertically.
 *  @param {HTMLElement[]} items
 *  @param {string} className - comma-separated CSS class names
 *  @returns {HTMLDivElement} */
export function divV(items: HTMLElement[], className?: string | null): HTMLDivElement

/** Div flex-box container that positions child elements horizontally.
 *  @param {HTMLElement[]} items
 *  @param {string} className - comma-separated CSS class names
 *  @returns {HTMLDivElement} */
export function divH(items: HTMLElement[], className?: string | null): HTMLDivElement

/** Renders content as a card. */
export function card(content: any): HTMLDivElement

export function loader(): any

export function setUpdateIndicator(element: HTMLElement, updating?: boolean): void

/**
 * Creates a button with the specified text, click handler, and tooltip
 * @param {string} text
 * @param {Function} handler
 * @param {string} tooltip
 * @returns {HTMLButtonElement}
 * */
export function button(text: string, handler: Function, tooltip?: string | null): HTMLButtonElement

export function bigButton(text: string, handler: Function, tooltip?: string | null): HTMLButtonElement

/**
 * Creates a combo popup with the specified icons and items
 * @param {HTMLElement} icon
 * @param {Array<string>} items
 * @param {Function} handler (item) => {...}
 * @returns {HTMLElement}
 * */
export function comboPopup(icon: HTMLElement, items: string[], handler: (item: any) => void): HTMLElement

/** Creates a visual table based on [map]. */
//TODO: what is this?
export function tableFromMap(map: any): any

/** Creates a visual element representing list of [items]. */
export function list(items: any[]): HTMLElement[]

/** Creates a [Dialog].
 * @returns {Dialog} */
export function dialog(title?: string | null): Dialog

/** Binds [item] with the [element]. It enables selecting it as a current object, drag-and-drop,
 * tooltip, and popup menu.
 * @returns Element. */
export function bind(item: any, element: Element): Element

/**
 * @param {string} value
 * @param {boolean} autoSize
 * @param {boolean} resizable
 * @returns Element */
export function textArea(value: string, autoSize?: boolean, resizable?: boolean): Element

/**
 * Creates a virtual list widget
 * @param {number} length - number of elements
 * @param {Function} renderer
 * @returns {VirtualView}
 */
export function virtualView(length: number, renderer: Function): VirtualView


export function popupMenu(items: any): void

export function tooltipHide(): void

//TODO: types...e and x variables?
//export function tooltip(e: HTMLElement, x: any): HTMLElement

export function tooltipShow(content: any, x: any, y: HTMLElement): HTMLElement

export function inputs(inputs: InputBase): HTMLElement;

/** Creates new nodes tree
 * @returns {TreeViewNode} */
export function tree(): TreeViewNode

export function intInput(name: string, value: number): InputBase

export function choiceInput<T>(name: string, selected: T, items: T[]): InputBase

export function multiChoiceInput<T>(name: string, value: T, items: T[]): InputBase

export function stringInput(name: string, value: string): InputBase

export function floatInput(name: string, value: number): InputBase

export function dateInput(name: string, value: number): InputBase

export function boolInput(name: string, value: boolean, callback: Function): InputBase

export function moleculeInput(name: string, value: any): InputBase

export function columnInput(name: string, table: any, value: any): InputBase

export function columnsInput(name: string, table: any): InputBase

/** UI Tools **/
export class tools {

    static handleResize(element: HTMLElement, onChanged: (width: number, height: number) => void): () => void;
}


/**
 * Override this class, and {@link register} an instance to integrate the platform with custom
 * types and objects.
 *
 * Samples: {@link https://public.datagrok.ai/js/samples/ui/meta/meta}
 * */
export class JsEntityMeta {

    /** Type of the object that this meta handles. */
    get type(): string

    static register(meta: JsEntityMeta): void

    /**
     * Override this method to check whether this meta class should handle the specified object.
     * @param x - specified object.
     * @returns {boolean}
     * */
    isApplicable(x: any): boolean

    /** String representation of the [item], by default item.toString().
     * @param x - item
     * @returns {string} */
    getCaption(x: any): string

    /** Renders icon for the item.
     * @param x - item
     * @returns {Element} */
    renderIcon(x: any): HTMLDivElement

    /** Renders markup for the item.
     * @param x - item
     * @returns {Element} */
    renderMarkup(x: any): HTMLDivElement

    /** Renders tooltip for the item.
     * @param x - item
     * @returns {Element} */
    renderTooltip(x: any): HTMLDivElement

    /** Renders card div for the item.
     * @param x - item
     * @returns {Element} */
    renderCard(x: any): HTMLDivElement

    /** Renders properties list for the item.
     * @param x - item
     * @returns {Element} */
    renderProperties(x: any): HTMLDivElement

    /** Renders view for the item.
     * @param x - item
     * @returns {Element} */
    renderView(x: any): HTMLDivElement

    /** Gets called once upon the registration of meta export class. */
    init(): void

    /**
     * Registers a function that takes applicable objects an the only argument.
     * It will be suggested to run in the context menu for that object, and
     * also in the "Actions" pane on the property panel.
     *
     * Samples: {@link https://public.datagrok.ai/js/samples/ui/docking/docking}
     *
     * @param {string} name - function name
     * @param run - a function that takes exactly one parameter
     * */
    registerParamFunc(name: string, run: (params: any) => void): void;
}
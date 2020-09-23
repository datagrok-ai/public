import {Observable} from "rxjs";
import {StreamSubscription} from "./events";
import {Property} from "./entities";
import {DataFrame} from "./dataframe";

/** Base class for controls that have a visual root and a set of properties. */
export class Widget {

    properties: Property[];
    
    /** @constructs Widget and initializes its root. */
    constructor()

    /** Widget's visual root.
     * @type {HTMLElement} */
    get root(): HTMLElement
    set root(r: HTMLElement);

    /** Creates a {@see Widget} from the specified React component. */
    static react(reactComponent: React.DOMElement<any, any> | Array<React.DOMElement<any, any>> | React.CElement<any, any> | Array<React.CElement<any, any>> | React.ReactElement | Array<React.ReactElement>): Widget
}

/** Base class for DataFrame-bound filtering controls */
export class Filter extends Widget {

    dataFrame: DataFrame | null;
    
    constructor()
}

export class DartWidget extends Widget {

    constructor(d: any)

    get root(): HTMLElement
}

/**
 * Accordion control with collapsible/expandable panes.
 * Samples: {@link https://public.datagrok.ai/js/samples/ui/accordion}
 * @extends {DartWidget}
 * */
export class Accordion extends DartWidget {

    /** @constructs Accordion */
    constructor(d: any)

    /** @type {AccordionPane[]} */
    get panes(): AccordionPane[]

    /** Creates a new instance of Accordion */
    static create(): Accordion

    /** Returns a pane with the specified name.
     * @param {string} name
     * @returns {AccordionPane} */
    getPane(name: string): AccordionPane

    addPane(name: string, getContent: Function, expanded?: boolean, before?: AccordionPane | null): void
}


/** A pane in the {@link Accordion} control. */
export class AccordionPane {
    constructor(d: any)

    /** Expanded state
     * @type {boolean} */
    get expanded(): boolean

    set expanded(v)

    /** @type {string} */
    get name(): string

    set name(name)
}


export class TabControl {
    constructor(d: any)

    get root(): HTMLDivElement

    get header(): HTMLDivElement

    get panes(): TabPane

    get currentPane(): TabPane

    set currentPane(v)

    static create(vertical?: boolean): TabControl

    getPane(name: string): TabPane

    addPane(name: string, getContent: Function, icon?: any | null): TabPane

    clear(): void

}


export class TabPane {
    constructor(d: any)

    get expanded(): boolean

    set expanded(v)

    get name(): string

    set name(name)
}


export class ToolboxPage {
    constructor(d: any)

    get accordion(): Accordion
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
export class Dialog {
    constructor(d: any)

    static create(title?: string): Dialog

    /**
     *  @param {Function} handler
     *  @returns {Dialog} */
    onOK(handler: () => void): Dialog

    /// Occurs when Dialog is closed
    onClose(callback: (v: any) => void): StreamSubscription

    /** @returns {Dialog} */
    show(): Dialog

    /** @returns {Dialog}
     * @param {boolean} fullScreen  */
    showModal(fullScreen: boolean): this


    /** Adds content to the dialog.
     * @param {HTMLElement | Widget | InputBase} content
     * @returns {Dialog} */
    add(content: HTMLElement | Widget | InputBase): Dialog
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

    constructor(d: any)

    static create(): Menu

    /** Creates a popup menu.
     * @returns {Menu} */
    static popup(): Menu

    /** Adds a menu group with the specified text.
     * @param {string} text*/
    group(text: string): Menu

    //TODO: string or?
    item(item: string, onClick: Function): Menu

    items(items: string[], onClick: Function): Menu

    /** Adds a separator line.
     *  @returns {Menu} */
    separator(): Menu

    /** Shows the menu.
     * @returns {Menu} */
    show(): Menu
}


/** Balloon-style visual notifications. */
export class Balloon {
    /** Shows information message (green background) */
    info(s: string): void

    /** Shows information message (red background) */
    error(s: string): void
}


/** Input control base. Could be used for editing {@link Property} values as well. */
export class InputBase {
    constructor(d: any, onChanged: (v: any) => void)

    get root(): HTMLElement

    get caption(): string

    get format(): string

    get captionLabel(): string

    //TODO: string or  HTML?
    get input(): string

    get nullable(): boolean

    set nullable(v)

    get value(): any

    set value(x)

    get stringValue(): string

    set stringValue(s)

    get readOnly(): boolean

    set readOnly(v)

    get enabled(): boolean

    set enabled(v)

    static forProperty(property: Property): InputBase

    /// Occurs when [value] is changed, either by user or programmatically.
    onChanged(callback: (v: any) => void): StreamSubscription

    /// Occurs when [value] is changed by user.
    onInput(callback: (v: any) => void): StreamSubscription

    save(): any

    load(s: any): any

    init(): any

    fireChanged(): any

    addCaption(caption: string): void

    addPatternMenu(pattern: any): void

    setTooltip(msg: string): any
}


export class ProgressIndicator {
    constructor(d: any)

    get percent(): number

    get description(): string

    set description(s)

    get onProgressUpdated(): Observable<any>

    get onLogUpdated(): Observable<any>

    static create(name?: string): ProgressIndicator

    update(percent: number, description: string): void
    
    log(line: any): void
}


export class TaskBarProgressIndicator extends ProgressIndicator {
    static create(name?: string): TaskBarProgressIndicator

    close(): any
}


export class TagEditor {
    constructor(d: any)

    get root(): HTMLElement

    get tags(): string

    set acceptsDragDrop(predicate: (...params: any[]) => boolean)

    set doDrop(action: Function)

    static create(): TagEditor

    //TODO: string or TagElement
    addTag(tag: TagElement, notify?: boolean): void

    removeTag(tag: TagElement): void

    clearTags(): void

    onChanged(callback: Function): StreamSubscription
}


export class TagElement {
    constructor(d: any)

    get tag(): string

    set tag(x)
}


type ColorType = number

/** Color-related routines. */
export class Color {

    static get categoricalPalette(): ColorType[]

    static get gray(): ColorType

    static get lightLightGray(): ColorType

    static get lightGray(): ColorType

    static get darkGray(): ColorType

    static get blue(): ColorType

    static get green(): ColorType

    static get darkGreen(): ColorType

    static get black(): ColorType

    static get yellow(): ColorType

    static get white(): ColorType

    static get red(): ColorType

    static get darkRed(): ColorType

    static get maroon(): ColorType

    static get olive(): ColorType

    static get orange(): ColorType

    static get darkOrange(): ColorType

    static get lightBlue(): ColorType

    static get darkBlue(): ColorType

    static get purple(): ColorType

    static get whitesmoke(): ColorType

    static get navy(): ColorType

    static get cyan(): ColorType

    static get filteredRows(): ColorType

    static get filteredOutRows(): ColorType

    static get selectedRows(): ColorType

    static get missingValueRows(): ColorType

    static get mouseOverRows(): ColorType

    static get currentRow(): ColorType

    static get histogramBar(): ColorType

    static get barChart(): ColorType

    static get scatterPlotMarker(): ColorType

    static get scatterPlotSelection(): ColorType

    static get scatterPlotZoom(): ColorType

    static get areaSelection(): ColorType

    static get rowSelection(): ColorType

    static get colSelection(): ColorType

    static get areaZoom(): ColorType

    static get gridWarningBackground(): ColorType

    static get success(): ColorType

    static get failure(): ColorType

    static r(c: number): number

    static g(c: number): number

    static b(c: number): number

    /** Returns i-th categorical color (looping over the palette if needed) */
    static getCategoricalColor(i: number): number

    /** Returns either black or white color, depending on which one would be most contrast to the specified [color] */
    static getContrastColor(color: number): ColorType

    static toRgb(color: ColorType): string

    static scale(x: number, min: number, max: number): number
}


/** Tree view node.
 * Sample: {@link https://public.datagrok.ai/js/samples/ui/tree-view}
 * */
export class TreeViewNode {
    /** @constructs {TreeView} */
    constructor(d: any)

    /** Visual root.
     * @type {HTMLElement} */
    get root(): HTMLElement

    /** Caption label.
     * @type {HTMLElement} */
    get captionLabel(): HTMLElement

    /** Check box.
     * @type {null|HTMLElement} */
    get checkBox(): HTMLElement | null

    /** Returns 'true' if checked
     * @returns {boolean} */
    get checked(): boolean

    set checked(checked)

    /** Returns node text.
     * @returns {string} */
    get text(): string

    /** Node value.
     * @type {Object}
     */
    get value(): Object

    set value(v)

    /** Gets all node items.
     * @returns {Array<TreeViewNode>} */
    get items(): TreeViewNode[]

    /** Creates new nodes tree
     * @returns {TreeViewNode} */
    static tree(): TreeViewNode

    /** Add new group to node.
     * @param {string} text
     * @param {Object} value
     * @param {boolean} expanded
     * @returns {TreeViewNode} */
    group(text: string, value?: Object | null, expanded?: boolean): TreeViewNode

    /** Add new item to node.
     * @param {string} text
     * @param {Object} value
     * @returns {TreeViewNode} */
    item(text: string, value?: Object | null): TreeViewNode

    /** Enables checkbox on node
     * @param {boolean} checked */
    enableCheckBox(checked?: boolean): void
}
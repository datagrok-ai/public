import {toDart, toJs} from "./wrappers";
import {_sub, observeStream} from "./events";

/** Base class for controls that have a visual root and a set of properties. */
export class Widget {

    /** @constructs Widget and initializes its root. */
    constructor() {
        /** @member {HTMLElement} */
        this._root = null;

        /** @member {Property[]}*/
        this.properties = [];
    }

    /** Widget's visual root.
     * @type {HTMLElement} */
    get root() { return this._root; }
    set root(r) { this._root = r; }



    /** Creates a {@see Widget} from the specified React component. */
    static react(reactComponent) {
        let widget = new Widget(ui.div());
        ReactDOM.render(reactComponent, widget.root);
        return widget;
    }
}


/** Base class for DataFrame-bound filtering controls */
export class Filter extends Widget {

    constructor() {
        super(ui.div());

        /** @member {DataFrame} */
        this.dataFrame = null;
    }
}


export class DartWidget extends Widget {

    constructor(d) {
        super();
        this.d = d;
    }

    get root() { return grok_TabControlBase_Get_Root(this.d); }
}

/**
 * Accordion control with collapsible/expandable panes.
 * Samples: {@link https://public.datagrok.ai/js/samples/ui/accordion}
 * @extends {DartWidget}
 * */
export class Accordion extends DartWidget {

    /** @constructs Accordion */
    constructor(d) { super(d); }

    /** Creates a new instance of Accordion */
    static create() { return toJs(grok_Accordion()); }

    /** @type {AccordionPane[]} */
    get panes() { return grok_TabControlBase_Get_Panes(this.d).map(toJs); }

    /** Returns a pane with the specified name.
     * @param {string} name
     * @returns {AccordionPane} */
    getPane(name) { return toJs(grok_TabControlBase_GetPane(this.d, name)); }

    addPane(name, getContent, expanded = false, before = null) {
        return toJs(grok_Accordion_AddPane(this.d, name, getContent, expanded, before !== null ? before.d : null));
    }
}


/** A pane in the {@link Accordion} control. */
export class AccordionPane {
    constructor(d) { this.d = d; }

    /** Expanded state
     * @type {boolean} */
    get expanded() { return grok_AccordionPane_Get_Expanded(this.d); }
    set expanded(v) { return grok_AccordionPane_Set_Expanded(this.d, v); }

    /** @type {string} */
    get name() { return grok_AccordionPane_Get_Name(this.d); }
    set name(name) { return grok_AccordionPane_Set_Name(this.d, name); }
}


export class TabControl {
    constructor(d) { this.d = d; }
    static create(vertical = false) { return toJs(grok_TabControl(vertical)); }

    get root() { return grok_TabControlBase_Get_Root(this.d); }
    get header() { return grok_TabControlBase_Get_Header(this.d); }
    get panes() { return grok_TabControlBase_Get_Panes(this.d).map(toJs); }
    getPane(name) { return toJs(grok_TabControlBase_GetPane(this.d, name)); }

    addPane(name, getContent, icon = null) {
        return toJs(grok_TabControlBase_AddPane(this.d, name, getContent, icon));
    }
    clear() { return toJs(grok_TabControlBase_Clear(this.d)); }

    get currentPane() { return grok_TabControlBase_Get_CurrentPane(this.d); }
    set currentPane(v) { return grok_TabControlBase_Set_CurrentPane(this.d, v.d); }

}


export class TabPane {
    constructor(d) { this.d = d; }

    get expanded() { return grok_AccordionPane_Get_Expanded(this.d); }
    set expanded(v) { return grok_AccordionPane_Set_Expanded(this.d, v); }

    get name() { return grok_AccordionPane_Get_Name(this.d); }
    set name(name) { return grok_AccordionPane_Set_Name(this.d, name); }
}


export class ToolboxPage {
    constructor(d) { this.d = d; }

    get accordion() { return toJs(grok_ToolboxPage_Get_Accordion(this.d)); }
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
    constructor(d) { this.d = d; }

    static create(title = '') { return new Dialog(grok_Dialog(title)); }

    /**
     *  @param {Function} handler
     *  @returns {Dialog} */
    onOK(handler) { grok_Dialog_OnOK(this.d, handler); return this; }

    /** @returns {Dialog} */
    show() { grok_Dialog_Show(this.d, false, false); return this; }

    /** @returns {Dialog}
     * @param {boolean} fullScreen  */
    showModal(fullScreen) { grok_Dialog_Show(this.d, true, fullScreen); return this; }

    /** Adds content to the dialog.
     * @param {HTMLElement | Widget | InputBase} content
     * @returns {Dialog} */
    add(content) { grok_Dialog_Add(this.d, toDart(ui.extract(content))); return this; }
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

    constructor(d) { this.d  = d; }

    static create() { return toJs(grok_Menu()); }

    /** Creates a popup menu.
     * @returns {Menu} */
    static popup() { return toJs(grok_Menu_Context()); }

    /** Adds a menu group with the specified text.
     * @param {string} text*/
    group(text) { return toJs(grok_Menu_Group(this.d, text)); }

    item(item, onClick) { return toJs(grok_Menu_Item(this.d, item, onClick)); }

    items(items, onClick) { return toJs(grok_Menu_Items(this.d, items, onClick)); }

    /** Adds a separator line.
     *  @returns {Menu} */
    separator() { return toJs(grok_Menu_Separator(this.d)); }

    /** Shows the menu.
     * @returns {Menu} */
    show() { return toJs(grok_Menu_Show(this.d)); }
}


/** Balloon-style visual notifications. */
export class Balloon {
    /** Shows information message (green background) */
    info(s) { grok_Balloon(s, 'info'); }

    /** Shows information message (red background) */
    error(s) { grok_Balloon(s, 'error'); }
}


/** Input control base. Could be used for editing {@link Property} values as well. */
export class InputBase {
    constructor(d, onChanged = null) {
        this.d = d;
        if (onChanged != null)
            this.onChanged((_) => onChanged(this.value));
    }

    get root() { return grok_InputBase_Get_Root(this.d); };
    get caption() { return grok_InputBase_Get_Caption(this.d); };
    get format() { return grok_InputBase_Get_Format(this.d); } ;
    get captionLabel() { return grok_InputBase_Get_CaptionLabel(this.d); };
    get input() { return grok_InputBase_Get_Input(this.d); };

    get nullable() { return grok_InputBase_Get_Nullable(this.d); };
    set nullable(v) { return grok_InputBase_Set_Nullable(this.d, v); };

    get value() { return grok_InputBase_Get_Value(this.d); };
    set value(x) { return grok_InputBase_Set_Value(this.d, x); };

    get stringValue() { return grok_InputBase_Get_StringValue(this.d); };
    set stringValue(s) { return grok_InputBase_Set_StringValue(this.d, s); };

    get readOnly() { return grok_InputBase_Get_ReadOnly(this.d); };
    set readOnly(v) { return grok_InputBase_Set_ReadOnly(this.d, v); };

    get enabled() { return grok_InputBase_Get_Enabled(this.d); };
    set enabled(v) { return grok_InputBase_Set_Enabled(this.d, v); };

    /// Occurs when [value] is changed, either by user or programmatically.
    onChanged(callback) { return _sub(grok_InputBase_OnChanged(this.d, callback)); }

    /// Occurs when [value] is changed by user.
    onInput(callback) { return _sub(grok_InputBase_OnInput(this.d, callback)); }

    save() { return grok_InputBase_Save(this.d); };
    load(s) { return grok_InputBase_Load(this.d, s); };

    init() { return grok_InputBase_Init(this.d); };
    fireChanged() { return grok_InputBase_FireChanged(this.d); };
    addCaption(caption) { grok_InputBase_AddCaption(this.d, caption); };
    addPatternMenu(pattern) { grok_InputBase_AddPatternMenu(this.d, pattern); }
    setTooltip(msg) { grok_InputBase_SetTooltip(this.d, msg); };

    static forProperty(property) { return toJs(grok_InputBase_ForProperty(property.d)); }
}


export class ProgressIndicator {
    constructor(d) { this.d = d; }

    static create() { return toJs(grok_ProgressIndicator_Create()); }

    get percent() { return grok_ProgressIndicator_Get_Percent(this.d); }

    get description() { return grok_ProgressIndicator_Get_Description(this.d); }
    set description(s) { grok_ProgressIndicator_Set_Description(this.d, s); }

    update(percent, description) { grok_ProgressIndicator_Update(this.d, percent, description); }
    log(line) { grok_ProgressIndicator_Log(this.d, line); }

    get onProgressUpdated() { return observeStream(grok_Progress_Updated(this.d)); }
    get onLogUpdated() { return observeStream(grok_Progress_Log_Updated(this.d)); }

}


export class TaskBarProgressIndicator extends ProgressIndicator {
    static create(name) { return toJs(grok_TaskBarProgressIndicator_Create(name)); }

    close() { return grok_TaskBarProgressIndicator_Close(this.d); }
}


export class TagEditor {
    constructor(d) { this.d = d; }

    static create() { return toJs(grok_TagEditor()); }

    get root() { return grok_TagEditor_Get_Root(this.d); }

    get tags() { return grok_TagEditor_Get_Tags(this.d); }

    addTag(tag, notify = true) { return grok_TagEditor_AddTag(this.d, tag, notify); }
    removeTag(tag) { grok_TagEditor_RemoveTag(this.d, tag); }
    clearTags() { grok_TagEditor_ClearTags(this.d); }

    set acceptsDragDrop(predicate) { grok_TagEditor_Set_AcceptsDragDrop(this.d, (x) => predicate(toJs(x, false))); };
    set doDrop(action) { grok_TagEditor_Set_DoDrop(this.d, (x) => action(toJs(x, false))); }

    onChanged(callback) { return _sub(grok_TagEditor_OnChanged(this.d, callback)); }
}


export class TagElement {
    constructor(d) { this.d = d; }

    get tag() { return grok_TagElement_Get_Tag(this.d); };
    set tag(x) { return grok_TagElement_Set_Tag(this.d, x); };
}


/** Color-related routines. */
export class Color {

    static r(c) { return (c >> 16) & 0xFF; }
    static g(c) { return (c >> 8) & 0xFF; }
    static b(c) { return c & 0xFF; }

    /** Returns i-th categorical color (looping over the palette if needed) */
    static getCategoricalColor(i) { return Color.categoricalPalette[i % Color.categoricalPalette.length]; }

    /** Returns either black or white color, depending on which one would be most contrast to the specified [color] */
    static getContrastColor(color) { return grok_Color_GetContrastColor(color); }

    static toRgb(color) { return color === null ? '': `rgb(${Color.r(color)},${Color.g(color)},${Color.b(color)})`; }

    static get categoricalPalette() { return grok_Color_CategoricalPalette(); }

    static scale(x, min, max) {
        return min === max ? min : (x - min) / (max - min);
    }

    static get gray() { return 0xFF808080; }
    static get lightLightGray() { return 0xFFF0F0F0; }
    static get lightGray() { return 0xFFD3D3D3; }
    static get darkGray() { return 0xFF838383; }
    static get blue() { return 0xFF0000FF; }
    static get green() { return 0xFF00FF00; }
    static get darkGreen() { return 0xFF006400; }
    static get black() { return 0xFF000000; }
    static get yellow() { return 0xFFFFFF00; }
    static get white() { return 0xFFFFFFFF; }
    static get red() { return 0xFFFF0000; }
    static get darkRed() { return 0xFF8b0000; }
    static get maroon() { return 0xFF800000; }
    static get olive() { return 0xFF808000; }
    static get orange() { return 0xFFFFA500; }
    static get darkOrange() { return 0xFFFF8C00; }
    static get lightBlue() { return 0xFFADD8E6; }
    static get darkBlue() { return 0xFF0000A0; }
    static get purple() { return 0xFF800080; }
    static get whitesmoke() { return 0xFFF5F5F5; }
    static get navy() { return 0xFF000080; }
    static get cyan() { return 0xFF00ffff; }

    static get filteredRows() { return 0xff1f77b4; }
    static get filteredOutRows() { return Color.lightLightGray; }
    static get selectedRows() { return Color.darkOrange; }
    static get missingValueRows() { return Color.filteredOutRows; }
    static get mouseOverRows() { return 0xFFAAAAAA; }
    static get currentRow() { return 0xFF38B738; }

    static get histogramBar() { return Color.filteredRows; }
    static get barChart() { return 0xFF24A221; }
    static get scatterPlotMarker() { return 0xFF40699c; }
    static get scatterPlotSelection() { return 0x80323232; }
    static get scatterPlotZoom() { return 0x80626200; }

    static get areaSelection() { return Color.lightBlue; }
    static get rowSelection() { return 0x60dcdca0; }
    static get colSelection() { return 0x60dcdca0; }
    static get areaZoom() { return 0x80323232; }

    static get gridWarningBackground() { return 0xFFFFB9A7; }

    static get success() { return 0xFF3cb173; }
    static get failure() { return 0xFFeb6767; }
}


/** Tree view node.
 * Sample: {@link https://public.datagrok.ai/js/samples/ui/tree-view}
 * */
export class TreeViewNode {
    /** @constructs {TreeView} */
    constructor(d) { this.d = d; }

    /** Creates new nodes tree
     * @returns {TreeViewNode} */
    static tree() { return toJs(grok_TreeViewNode_Tree()); }

    /** Visual root.
     * @type {HTMLElement} */
    get root() { return grok_TreeViewNode_Root(this.d); }

    /** Caption label.
     * @type {HTMLElement} */
    get captionLabel() { grok_TreeViewNode_CaptionLabel(this.d); }

    /** Check box.
     * @type {null|HTMLElement} */
    get checkBox() { return grok_TreeViewNode_CheckBox(this.d); }

    /** Returns 'true' if checked
     * @returns {boolean} */
    get checked() { return grok_TreeViewNode_Get_Checked(this.d); }
    set checked(checked) { return grok_TreeViewNode_Set_Checked(this.d, checked); }

    /** Returns node text.
     * @returns {string} */
    get text() { return grok_TreeViewNode_Text(this.d); }

    /** Node value.
     * @type {Object}
     */
    get value() { return grok_TreeViewNode_Get_Value(this.d) };
    set value(v) { grok_TreeViewNode_Set_Value(this.d, v) };

    /** Gets all node items.
     * @returns {Array<TreeViewNode>} */
    get items() { return grok_TreeViewNode_Items(this.d).map(i => toJs(i)); }

    /** Add new group to node.
     * @param {string} text
     * @param {Object} value
     * @param {boolean} expanded
     * @returns {TreeViewNode} */
    group(text, value = null, expanded = true) { return toJs(grok_TreeViewNode_Group(this.d, text, value, expanded)); }

    /** Add new item to node.
     * @param {string} text
     * @param {Object} value
     * @returns {TreeViewNode} */
    item(text, value = null) { return toJs(grok_TreeViewNode_Item(this.d, text, value)); }

    /** Enables checkbox on node
     * @param {boolean} checked */
    enableCheckBox(checked = false) { grok_TreeViewNode_EnableCheckBox(this.d, checked); }
}
